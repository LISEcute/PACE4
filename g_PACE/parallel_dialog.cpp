#include "parallel_dialog.h"

#include <QCloseEvent>
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDateTime>
#include <QDesktopServices>
#include <QDir>
#include <QDirIterator>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontMetrics>
#include <QHeaderView>
#include <QHBoxLayout>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QProcess>
#include <QProcessEnvironment>
#include <QProgressBar>
#include <QPushButton>
#include <QRegularExpression>
#include <QSaveFile>
#include <QSet>
#include <QSpinBox>
#include <QStandardPaths>
#include <QTableWidget>
#include <QThread>
#include <QTextBrowser>
#include <QTimer>
#include <QUrl>
#include <QVBoxLayout>
#include <algorithm>
#include <cmath>
#ifdef Q_OS_WIN
#include <qt_windows.h>
#endif

namespace {
constexpr quint32 modulus = 2147483647u;
constexpr quint32 mask = 123459876u;
quint32 effectiveSeed(quint32 seed) {
    const quint32 state = seed ^ mask;
    return state == 0 || state == modulus ? 1 : state;
}
QString hash(const QByteArray& bytes) {
    return QString::fromLatin1(QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
}
QString now() { return QDateTime::currentDateTimeUtc().toString(Qt::ISODateWithMs); }
bool readBytes(const QString& path, QByteArray& bytes, QString& error) {
    QFile file(path);
    if(!file.open(QIODevice::ReadOnly)) { error = path + ": " + file.errorString(); return false; }
    bytes = file.readAll();
    if(file.error() != QFileDevice::NoError) { error = path + ": " + file.errorString(); return false; }
    return true;
}
bool writeBytes(const QString& path, const QByteArray& bytes, QString& error) {
    QSaveFile file(path);
    if(!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() || !file.commit()) {
        error = path + ": " + file.errorString(); return false;
    }
    return true;
}
bool headerCount(const QByteArray& input, int& count) {
    QByteArray header = input.left(input.indexOf('\n') < 0 ? input.size() : input.indexOf('\n'));
    if(header.startsWith("\xEF\xBB\xBF")) header.remove(0, 3);
    const QStringList fields = QString::fromUtf8(header).trimmed().split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
    if(fields.size() != 10 && fields.size() != 11) return false;
    if(!QRegularExpression("^[0-9]+$").match(fields[0]).hasMatch()) return false;
    bool ok = false;
    count = fields[0].toInt(&ok);
    return ok && count >= 1 && count <= 1000000;
}
QTableWidgetItem* textItem(const QString& text) {
    auto* item = new QTableWidgetItem(text);
    item->setFlags(item->flags() & ~Qt::ItemIsEditable);
    return item;
}
}

ParallelDialog::ParallelDialog(QWidget* parent, const QString& executable)
    : QDialog(parent), m_executable(executable.isEmpty() ? QCoreApplication::applicationFilePath() : QFileInfo(executable).absoluteFilePath())
{
    setWindowTitle(tr("Parallel calculations"));
    setObjectName("parallelDialog");
    setModal(true);
    resize(940, 640);
    auto* layout = new QVBoxLayout(this);
    auto* intro = new QLabel(tr("Add saved .pace files. Duplicate an input to run another replica with a different seed. "
                               "Save changes in the main window first. Each calculation keeps its original cascade count."), this);
    intro->setWordWrap(true);
    layout->addWidget(intro);
    auto* inputButtons = new QHBoxLayout;
    m_add = new QPushButton(tr("Add inputs…"), this); m_add->setObjectName("addInputs");
    m_duplicate = new QPushButton(tr("Duplicate selected"), this); m_duplicate->setObjectName("duplicateInputs");
    m_remove = new QPushButton(tr("Remove selected"), this); m_remove->setObjectName("removeInputs");
    inputButtons->addWidget(m_add); inputButtons->addWidget(m_duplicate); inputButtons->addWidget(m_remove); inputButtons->addStretch();
    layout->addLayout(inputButtons);
    m_table = new QTableWidget(0, 4, this); m_table->setObjectName("inputTable");
    m_table->setHorizontalHeaderLabels({tr("Input file"), tr("Cascades"), tr("Seed"), tr("Status")});
    m_table->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_table->setSelectionMode(QAbstractItemView::ExtendedSelection);
    m_table->setSortingEnabled(false);
    m_table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    for(int column = 1; column < 4; ++column) m_table->horizontalHeader()->setSectionResizeMode(column, QHeaderView::ResizeToContents);
    m_table->horizontalHeader()->setSectionResizeMode(2, QHeaderView::Fixed);
    m_table->setColumnWidth(2, std::max(170, fontMetrics().horizontalAdvance("2147483646") + 70));
    layout->addWidget(m_table, 2);
    auto* controls = new QHBoxLayout;
    m_workers = new QSpinBox(this); m_workers->setObjectName("workerCount"); m_workers->setRange(1, 64);
    m_workers->setValue(std::min(2, std::max(1, QThread::idealThreadCount())));
    m_workers->setMinimumWidth(std::max(90, m_workers->fontMetrics().horizontalAdvance("64") + 70));
    m_workers->setToolTip(tr("Number of independent calculations running at once. Choose 1 for a serial run."));
    controls->addWidget(new QLabel(tr("Workers:"), this)); controls->addWidget(m_workers);
    m_timeout = new QSpinBox(this); m_timeout->setObjectName("timeoutSeconds"); m_timeout->setRange(0, 86400);
    m_timeout->setSpecialValueText(tr("No limit")); m_timeout->setSuffix(tr(" s"));
    m_timeout->setMinimumWidth(std::max(150, m_timeout->fontMetrics().horizontalAdvance(tr("No limit")) + 70));
    controls->addSpacing(20); controls->addWidget(new QLabel(tr("Timeout per job:"), this)); controls->addWidget(m_timeout); controls->addStretch();
    layout->addLayout(controls);
    auto* outputRow = new QHBoxLayout;
    m_outputEdit = new QLineEdit(proposedOutput(QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation) + "/PACE4 runs"), this);
    m_outputEdit->setObjectName("outputDirectory"); m_outputEdit->setPlaceholderText(tr("New results folder"));
    m_browse = new QPushButton(tr("Browse…"), this); m_browse->setObjectName("browseOutput");
    outputRow->addWidget(new QLabel(tr("Results folder:"), this)); outputRow->addWidget(m_outputEdit, 1); outputRow->addWidget(m_browse);
    layout->addLayout(outputRow);
    auto* outputHint = new QLabel(tr("Choose a new folder for this batch. Existing folders are kept unchanged."), this);
    layout->addWidget(outputHint);
    m_progress = new QProgressBar(this); m_progress->setObjectName("progressBar"); m_progress->setRange(0, 1); m_progress->setValue(0);
    layout->addWidget(m_progress);
    m_status = new QLabel(tr("Ready. Add one or more input files."), this); m_status->setObjectName("statusLabel"); m_status->setWordWrap(true);
    layout->addWidget(m_status);
    m_log = new QPlainTextEdit(this); m_log->setObjectName("logOutput"); m_log->setReadOnly(true); m_log->setMaximumBlockCount(2000);
    layout->addWidget(m_log, 1);
    auto* actions = new QHBoxLayout;
    m_open = new QPushButton(tr("Open results"), this); m_open->setObjectName("openResults"); m_open->setEnabled(false);
    m_start = new QPushButton(tr("Start batch"), this); m_start->setObjectName("startBatch"); m_start->setDefault(true);
    m_cancel = new QPushButton(tr("Cancel batch"), this); m_cancel->setObjectName("cancelBatch"); m_cancel->setEnabled(false);
    m_close = new QPushButton(tr("Close"), this); m_close->setObjectName("closeDialog");
    auto* helpButton = new QPushButton(tr("Help"), this);
    helpButton->setObjectName("parallelHelp");
    helpButton->setAutoDefault(false);
    helpButton->setShortcut(QKeySequence(Qt::Key_F1));
    helpButton->setToolTip(tr("How to run independent calculations (F1)"));
    auto* helpDialog = new QDialog(this);
    helpDialog->setObjectName("parallelHelpDialog");
    helpDialog->setWindowTitle(tr("Parallel calculations help"));
    helpDialog->resize(720, 620);
    auto* helpLayout = new QVBoxLayout(helpDialog);
    auto* helpText = new QTextBrowser(helpDialog);
    helpText->setObjectName("parallelHelpText");
    helpText->setHtml(tr(
        "<h2>Using the native parallel interface</h2>"
        "<p>Open <b>File &gt; Parallel calculations...</b>, click <b>Parallel...</b> on the toolbar, "
        "or press <b>Ctrl+Shift+P</b>. The dialog runs complete, independent calculations through "
        "Qt child processes. Python is not required. Save edits in the main window first: the queue "
        "reads saved files, and an existing current input is added automatically.</p>"
        "<ol><li><p>Click <b>Add inputs...</b> and select one or more saved .pace files. Use "
        "<b>Duplicate selected</b> to create another replica, or <b>Remove selected</b> to remove a row.</p></li>"
        "<li><p>Review the cascade counts and editable <b>seeds</b>. Each row needs a distinct effective "
        "RNG starting state. Keep its input and seed unchanged to reproduce its scientific output.</p></li>"
        "<li><p>Set <b>Workers</b> to the maximum number of simultaneous calculations. The range is 1 to 64; "
        "the initial value is 2 when at least two logical processors are available. One worker provides "
        "a serial reference. Higher counts require more memory.</p></li>"
        "<li><p>Choose a new <b>Results folder</b>. <b>Browse</b> selects a parent and proposes a new batch "
        "subfolder. Existing folders are refused. Optionally set <b>Timeout per job</b> in seconds; "
        "<b>No limit</b> disables the timeout.</p></li>"
        "<li><p>Click <b>Start batch</b>. Watch per-job status, completed-job progress and the log. "
        "<b>Cancel batch</b> stops active children and cancels queued jobs; closing the calculations "
        "dialog during a run also cancels and waits for cleanup.</p></li>"
        "<li><p>Select a completed row and click <b>Open results</b> to open its HTML report. With no row "
        "selected, this opens the batch folder. <b>Run again</b> retains row identities and seeds while "
        "proposing a fresh folder; input files are read again, so preserve them for replay.</p></li></ol>"
        "<p><b>Each job keeps its complete cascade count and probability cutoff.</b> Workers controls "
        "how many separate jobs run at once; it does not split one calculation. Changing the worker "
        "count does not change a job's seed. Different physical inputs must not be pooled as replicas.</p>"
        "<p>You can leave this help open during a batch. Closing this help window does not cancel calculations.</p>"));
    helpLayout->addWidget(helpText);
    auto* helpClose = new QDialogButtonBox(QDialogButtonBox::Close, helpDialog);
    helpLayout->addWidget(helpClose);
    connect(helpClose, &QDialogButtonBox::rejected, helpDialog, &QDialog::close);
    connect(helpButton, &QPushButton::clicked, helpDialog, [helpDialog]() {
        helpDialog->show(); helpDialog->raise(); helpDialog->activateWindow();
    });
    actions->addWidget(helpButton);
    actions->addWidget(m_open); actions->addStretch(); actions->addWidget(m_start); actions->addWidget(m_cancel); actions->addWidget(m_close);
    layout->addLayout(actions);
    m_timer = new QTimer(this); m_timer->setInterval(100);
    connect(m_timer, &QTimer::timeout, this, &ParallelDialog::checkTimeouts);
    connect(m_add, &QPushButton::clicked, this, [this]() {
        addInputFiles(QFileDialog::getOpenFileNames(this, tr("Add PACE4 inputs"), QString(), tr("PACE4 inputs (*.pace);;All files (*)")));
    });
    connect(m_duplicate, &QPushButton::clicked, this, [this]() {
        QStringList files;
        for(const auto& row : m_table->selectionModel()->selectedRows()) files << m_jobs[row.row()].source;
        addInputFiles(files);
    });
    connect(m_remove, &QPushButton::clicked, this, [this]() {
        QList<int> rows;
        for(const auto& row : m_table->selectionModel()->selectedRows()) rows << row.row();
        std::sort(rows.begin(), rows.end(), std::greater<int>());
        for(int row : rows) m_jobs.removeAt(row);
        refreshRows();
    });
    connect(m_browse, &QPushButton::clicked, this, [this]() {
        const QString parentFolder = QFileDialog::getExistingDirectory(this, tr("Choose parent folder for results"), QFileInfo(m_outputEdit->text()).absolutePath());
        if(!parentFolder.isEmpty()) m_outputEdit->setText(proposedOutput(parentFolder));
    });
    connect(m_start, &QPushButton::clicked, this, [this]() { startBatch(); });
    connect(m_cancel, &QPushButton::clicked, this, &ParallelDialog::cancelBatch);
    connect(m_close, &QPushButton::clicked, this, &ParallelDialog::reject);
    connect(m_open, &QPushButton::clicked, this, &ParallelDialog::openSelectedResults);
    connect(m_table, &QTableWidget::cellDoubleClicked, this, [this](int row, int column) {
        if(column == 3 && !m_jobs[row].folder.isEmpty()) openSelectedResults();
    });
}

ParallelDialog::~ParallelDialog() {
    // Application shutdown is the only path that can bypass the asynchronous close guard.
    for(auto& job : m_jobs) if(job.process) { job.process->disconnect(this); job.process->kill(); }
    for(auto& job : m_jobs) if(job.process) job.process->waitForFinished(1000);
}
QString ParallelDialog::proposedOutput(const QString& parent) const {
    return QDir(parent).filePath("run-" + QDateTime::currentDateTime().toString("yyyyMMdd-hhmmss-zzz"));
}
void ParallelDialog::addInputFiles(const QStringList& paths) {
    if(m_running) return;
    QSet<quint32> used;
    for(const auto& job : m_jobs) used.insert(effectiveSeed(job.seed));
    for(const auto& path : paths) {
        while(used.contains(effectiveSeed(m_nextSeed))) m_nextSeed = (quint64(m_nextSeed) + 104729u) % (modulus - 1u) + 1u;
        Job job;
        job.id = QString("job-%1").arg(m_nextIdentity++, 4, 10, QChar('0'));
        job.source = QFileInfo(path).absoluteFilePath(); job.seed = m_nextSeed;
        used.insert(effectiveSeed(job.seed));
        m_nextSeed = (quint64(m_nextSeed) + 104729u) % (modulus - 1u) + 1u;
        m_jobs.append(job);
    }
    refreshRows();
}
void ParallelDialog::setOutputDirectory(const QString& path) { if(!m_running) m_outputEdit->setText(path); }
void ParallelDialog::setWorkerCount(int count) { if(!m_running) m_workers->setValue(count); }
void ParallelDialog::setTimeoutSeconds(int seconds) { if(!m_running) m_timeout->setValue(seconds); }
void ParallelDialog::setJobSeed(int row, quint32 seed) {
    if(m_running || row < 0 || row >= m_jobs.size()) return;
    m_jobs[row].seed = seed; refreshRows();
}
void ParallelDialog::refreshRows() {
    m_table->setRowCount(0);
    for(int row = 0; row < m_jobs.size(); ++row) {
        Job& job = m_jobs[row];
        m_table->insertRow(row);
        auto* file = textItem(QFileInfo(job.source).fileName()); file->setToolTip(job.source + "\n" + job.id);
        m_table->setItem(row, 0, file);
        QByteArray bytes; QString ignored; int count = 0;
        const bool valid = readBytes(job.source, bytes, ignored) && headerCount(bytes, count);
        m_table->setItem(row, 1, textItem(valid ? QString::number(count) : tr("Invalid")));
        auto* seed = new QSpinBox(m_table); seed->setObjectName(QString("jobSeed_%1").arg(row)); seed->setRange(1, int(modulus - 1)); seed->setValue(int(job.seed));
        seed->setMinimumWidth(std::max(170, seed->fontMetrics().horizontalAdvance("2147483646") + 70));
        seed->setToolTip(tr("Keep this seed to reproduce this calculation. Every row needs a different RNG starting state."));
        const QString id = job.id;
        connect(seed, &QSpinBox::valueChanged, this, [this, id](int value) {
            for(auto& item : m_jobs) if(item.id == id) { item.seed = quint32(value); break; }
        });
        m_table->setCellWidget(row, 2, seed); m_table->setItem(row, 3, textItem(job.state));
    }
}
void ParallelDialog::setEditingEnabled(bool enabled) {
    for(auto* button : {m_add, m_duplicate, m_remove, m_browse, m_start}) button->setEnabled(enabled);
    m_workers->setEnabled(enabled); m_timeout->setEnabled(enabled); m_outputEdit->setEnabled(enabled);
    for(int row = 0; row < m_table->rowCount(); ++row) m_table->cellWidget(row, 2)->setEnabled(enabled);
    m_cancel->setEnabled(!enabled); m_close->setText(enabled ? tr("Close") : tr("Cancel and close"));
}
void ParallelDialog::appendLog(const QString& message) { m_log->appendPlainText(message); }
void ParallelDialog::showError(const QString& message) {
    m_lastError = message; m_status->setText(message); appendLog(message); emit validationFailed(message);
}
void ParallelDialog::setState(int row, const QString& state) {
    m_jobs[row].state = state; m_jobs[row].record["status"] = state;
    m_table->item(row, 3)->setText(state); emit jobStateChanged(row, state);
}

bool ParallelDialog::startBatch() {
    if(m_running) return false;
    m_lastError.clear();
    if(m_jobs.isEmpty()) { showError(tr("Add at least one saved .pace input.")); return false; }
    if(m_outputEdit->text().trimmed().isEmpty()) { showError(tr("Choose a new results folder.")); return false; }
    const QString output = QDir::cleanPath(QFileInfo(m_outputEdit->text().trimmed()).absoluteFilePath());
    if(QFileInfo::exists(output)) { showError(tr("The results folder already exists. Choose a new folder.")); return false; }
    QString error;
    if(!readBytes(m_executable, m_binary, error) || !readBytes(QFileInfo(m_executable).absolutePath() + "/lisecfg/AME_DB.sqlite", m_database, error)) {
        showError(error); return false;
    }
    QSet<quint32> used;
    for(auto& job : m_jobs) {
        if(job.seed < 1 || job.seed >= modulus || used.contains(effectiveSeed(job.seed))) {
            showError(tr("Every job needs a unique seed from 1 to 2147483646. Two seeds may also refer to the same masked RNG starting state.")); return false;
        }
        used.insert(effectiveSeed(job.seed));
        if(!readBytes(job.source, job.input, error)) { showError(error); return false; }
        if(!headerCount(job.input, job.ncasc)) { showError(tr("%1: expected a PACE header with 10 or 11 fields and 1–1000000 cascades.").arg(job.source)); return false; }
    }
    const QFileInfo outputInfo(output);
    QDir parent(outputInfo.absolutePath());
    if(!QDir().mkpath(parent.absolutePath()) || !parent.mkdir(outputInfo.fileName())) {
        showError(tr("Cannot create the new results folder: %1").arg(output)); return false;
    }
    m_output = output; m_running = true; m_stopping = m_failed = m_closeWhenFinished = m_persistFailed = false;
    m_nextJob = m_active = 0; m_totalClock.start(); m_executionClock.invalidate();
    m_log->clear(); setEditingEnabled(false); m_open->setEnabled(true); m_progress->setRange(0, m_jobs.size()); m_progress->setValue(0);
    m_summary = QJsonObject{{"version", 1}, {"status", "preparing"}, {"started_utc", now()},
        {"executable", m_executable}, {"executable_sha256", hash(m_binary)}, {"database_sha256", hash(m_database)},
        {"workers", m_workers->value()}, {"timeout_seconds", m_timeout->value()},
        {"seed_policy", "Explicit per-job seed; independent of worker count."},
        {"cutoff_policy", "Original NCASC preserved byte-for-byte; no cascade splitting."}};
    QJsonArray manifestJobs;
    for(int row = 0; row < m_jobs.size(); ++row) {
        Job& job = m_jobs[row]; job.folder = m_output + "/jobs/" + job.id;
        job.record = QJsonObject{{"id", job.id}, {"input_source", job.source}, {"input_sha256", hash(job.input)},
            {"seed", qint64(job.seed)}, {"effective_rng_initial_state", qint64(effectiveSeed(job.seed))},
            {"ncasc", job.ncasc}, {"effective_discriminator", 1.0 / std::max(500.0, 3.0 * std::sqrt(double(job.ncasc)) + 10.0)},
            {"folder", job.folder}, {"seconds", 0.0}, {"outputs", QJsonObject()}};
        setState(row, "pending");
        manifestJobs.append(QJsonObject{{"id", job.id}, {"input", job.source}, {"seed", qint64(job.seed)}});
    }
    if(!writeSummary() || !writeBytes(m_output + "/manifest.json", QJsonDocument(QJsonObject{{"version", 1}, {"jobs", manifestJobs}}).toJson(), error)) {
        stopBatch(error.isEmpty() ? m_lastError : error, true); return false;
    }
    // Snapshotting is completed before any child starts. All byte copies are private.
    for(auto& job : m_jobs) {
        const QString runtime = job.folder + "/runtime";
        const QString localExe = runtime + "/" + QFileInfo(m_executable).fileName();
        if(!QDir().mkpath(runtime + "/lisecfg") || !QDir().mkpath(runtime + "/files") ||
           !writeBytes(localExe, m_binary, error) || !writeBytes(runtime + "/lisecfg/AME_DB.sqlite", m_database, error) ||
           !writeBytes(runtime + "/lisepp.ini", "[Version]\nWorkInROOT=1\n[font]\nsize=9\nscaling=0\n", error) ||
           !writeBytes(job.folder + "/input.pace", job.input, error)) {
            stopBatch(error.isEmpty() ? tr("Could not prepare private runtime: %1").arg(runtime) : error, true); return false;
        }
        QFile::setPermissions(localExe, QFile::permissions(m_executable));
        job.record["runtime_executable"] = localExe;
    }
    m_summary["setup_seconds"] = m_totalClock.elapsed() / 1000.0;
    m_summary["status"] = "running"; m_executionClock.start();
    if(!writeSummary()) { stopBatch(m_lastError, true); return false; }
    appendLog(tr("Starting %1 independent calculations with %2 workers.").arg(m_jobs.size()).arg(m_workers->value()));
    m_timer->start(); QTimer::singleShot(0, this, &ParallelDialog::schedule);
    return true;
}

void ParallelDialog::schedule() {
    if(!m_running || m_stopping) { completeIfStopped(); return; }
    while(m_active < m_workers->value() && m_nextJob < m_jobs.size() && !m_stopping) launch(m_nextJob++);
    checkTimeouts(); completeIfStopped();
}
void ParallelDialog::launch(int row) {
    Job& job = m_jobs[row];
    auto* process = new QProcess(this); job.process = process;
    const QString originalDir = QFileInfo(m_executable).absolutePath();
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    for(const QString& key : env.keys()) if(key.startsWith("PACE_BENCHMARK_")) env.remove(key);
    env.insert("PACE4_SEED", QString::number(job.seed));
#ifdef Q_OS_WIN
    env.insert("QT_QPA_PLATFORM", QFileInfo::exists(originalDir + "/platforms/qoffscreen.dll") ? "offscreen" : "windows");
#else
    env.insert("QT_QPA_PLATFORM", "offscreen");
#endif
    const QString separator(QDir::listSeparator());
    env.insert("PATH", originalDir + separator + env.value("PATH"));
    env.insert("QT_PLUGIN_PATH", originalDir + separator + originalDir + "/plugins" + separator + env.value("QT_PLUGIN_PATH"));
    if(QFileInfo(originalDir + "/platforms").isDir()) env.insert("QT_QPA_PLATFORM_PLUGIN_PATH", originalDir + "/platforms");
    process->setProcessEnvironment(env); process->setWorkingDirectory(job.folder);
    process->setStandardInputFile(QProcess::nullDevice());
    process->setStandardOutputFile(job.folder + "/stdout.log"); process->setStandardErrorFile(job.folder + "/stderr.log");
    process->setProgram(job.record["runtime_executable"].toString()); process->setArguments({"input.pace", "-s"});
#ifdef Q_OS_WIN
    process->setCreateProcessArgumentsModifier([](QProcess::CreateProcessArguments* args) { args->flags |= CREATE_NO_WINDOW; });
#endif
    connect(process, &QProcess::started, this, [this, row]() {
        m_jobs[row].record["pid"] = qint64(m_jobs[row].process->processId());
        if(!writeSummary()) stopBatch(m_lastError, true);
    });
    connect(process, &QProcess::errorOccurred, this, [this, row](QProcess::ProcessError error) {
        if(error == QProcess::FailedToStart && m_jobs[row].process) {
            m_jobs[row].record["error"] = m_jobs[row].process->errorString(); processFinished(row, -1, true);
        }
    });
    connect(process, qOverload<int, QProcess::ExitStatus>(&QProcess::finished), this, [this, row](int code, QProcess::ExitStatus status) {
        processFinished(row, code, status == QProcess::CrashExit);
    });
    ++m_active; job.elapsed.start(); job.record["started_utc"] = now();
    job.record["command"] = QJsonArray{process->program(), "input.pace", "-s"}; setState(row, "running");
    appendLog(tr("%1 started · seed %2 · %3 cascades").arg(job.id).arg(job.seed).arg(job.ncasc));
#ifdef Q_OS_WIN
    const UINT previousMode = SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);
#endif
    process->start();
#ifdef Q_OS_WIN
    SetErrorMode(previousMode);
#endif
    if(!writeSummary()) stopBatch(m_lastError, true);
}
void ParallelDialog::collectOutputs(Job& job) {
    QJsonObject outputs;
    QStringList paths;
    QDir folder(job.folder);
    for(const QString& file : folder.entryList({"input.*", "stdout.log", "stderr.log"}, QDir::Files))
        if(file != "input.pace") paths << folder.filePath(file);
    QDirIterator trace(job.folder + "/runtime/files", QDir::Files, QDirIterator::Subdirectories);
    while(trace.hasNext()) paths << trace.next();
    for(const QString& path : paths) {
        QFile file(path);
        if(!file.open(QIODevice::ReadOnly)) continue;
        QCryptographicHash digest(QCryptographicHash::Sha256); digest.addData(&file);
        outputs[folder.relativeFilePath(path)] = QJsonObject{{"sha256", QString::fromLatin1(digest.result().toHex())}, {"bytes", file.size()}};
    }
    job.record["outputs"] = outputs;
}
void ParallelDialog::processFinished(int row, int exitCode, bool crashed) {
    Job& job = m_jobs[row];
    if(!job.process) return;
    job.process->disconnect(this); job.process->deleteLater(); job.process = nullptr; --m_active;
    job.record["returncode"] = exitCode; job.record["finished_utc"] = now(); job.record["seconds"] = job.elapsed.elapsed() / 1000.0;
    QStringList missing;
    for(const QString& name : {QString("input.cs4"), QString("input.html"), QString("input.rtf")}) {
        const QFileInfo output(job.folder + "/" + name);
        if(!output.isFile() || output.size() == 0) missing << name;
    }
    const bool failed = crashed || exitCode != 0 || !missing.isEmpty();
    if(job.state == "timed_out") setState(row, "timed_out");
    else if(m_stopping) setState(row, "canceled");
    else if(failed) {
        setState(row, "failed");
        if(!job.record.contains("error")) job.record["error"] = tr("Exit %1; missing or empty outputs: %2").arg(exitCode).arg(missing.join(", "));
    } else setState(row, "succeeded");
    collectOutputs(job);
    appendLog(tr("%1: %2 (%3 s)").arg(job.id, job.state).arg(job.elapsed.elapsed() / 1000.0, 0, 'f', 2));
    if(job.record.contains("error")) appendLog(job.record["error"].toString());
    if(failed && !m_stopping) {
        QFile errors(job.folder + "/stderr.log");
        if(errors.open(QIODevice::ReadOnly)) appendLog(QString::fromLocal8Bit(errors.read(4000)));
        stopBatch(tr("%1 failed. Other jobs are being stopped; logs and partial results are kept.").arg(job.id), true);
    }
    if(!writeSummary() && !m_stopping) stopBatch(m_lastError, true);
    QTimer::singleShot(0, this, &ParallelDialog::schedule);
}
void ParallelDialog::cancelBatch() {
    if(m_running) stopBatch(tr("Batch canceled. Partial results and logs are kept."), false);
}
void ParallelDialog::stopBatch(const QString& reason, bool failure) {
    if(!m_running) return;
    if(failure) m_failed = true;
    if(!m_stopping) {
        m_stopping = true; m_lastError = failure ? reason : QString(); m_summary["error"] = reason; appendLog(reason); m_status->setText(reason);
        for(int row = 0; row < m_jobs.size(); ++row) {
            Job& job = m_jobs[row];
            if(job.state == "pending") { job.record["error"] = reason; setState(row, "canceled"); }
            if(job.process) {
                const QPointer<QProcess> process = job.process;
                process->terminate();
                QTimer::singleShot(1500, this, [process]() { if(process && process->state() != QProcess::NotRunning) process->kill(); });
            }
        }
    }
    writeSummary(); completeIfStopped();
}
void ParallelDialog::checkTimeouts() {
    if(!m_running) return;
    int completed = 0;
    for(int row = 0; row < m_jobs.size(); ++row) {
        Job& job = m_jobs[row];
        if(job.state != "pending" && job.state != "running" && !job.process) ++completed;
        if(!m_stopping && job.process && m_timeout->value() > 0 && job.elapsed.elapsed() >= qint64(m_timeout->value()) * 1000) {
            job.record["error"] = tr("Job exceeded the %1 second timeout.").arg(m_timeout->value()); setState(row, "timed_out");
            stopBatch(job.record["error"].toString(), true);
        }
    }
    m_progress->setValue(completed);
    if(!m_stopping) m_status->setText(tr("%1 of %2 complete · %3 running").arg(completed).arg(m_jobs.size()).arg(m_active));
}
void ParallelDialog::completeIfStopped() {
    if(m_running && m_active == 0 && (m_stopping || m_nextJob >= m_jobs.size())) finishBatch();
}
void ParallelDialog::finishBatch() {
    if(!m_running) return;
    m_timer->stop(); m_running = false;
    const QString status = m_failed ? "failed" : (m_stopping ? "canceled" : "succeeded");
    m_summary["status"] = status; m_summary["finished_utc"] = now(); m_summary["total_seconds"] = m_totalClock.elapsed() / 1000.0;
    m_summary["execution_seconds"] = m_executionClock.isValid() ? m_executionClock.elapsed() / 1000.0 : 0.0;
    if(!m_summary.contains("setup_seconds")) m_summary["setup_seconds"] = m_totalClock.elapsed() / 1000.0;
    if(!writeSummary()) { m_failed = true; m_summary["status"] = "failed"; }
    m_progress->setValue(m_jobs.size()); setEditingEnabled(true);
    m_status->setText(m_failed ? tr("Batch failed. See the job status and log below.") : m_stopping ? tr("Batch canceled. Partial results are available.") : tr("All calculations completed. Select a row and open its results."));
    m_start->setText(tr("Run again"));
    // A repeat keeps identities/seeds while proposing a fresh output directory.
    m_outputEdit->setText(proposedOutput(QFileInfo(m_output).absolutePath()));
    m_binary.clear(); m_database.clear(); for(auto& job : m_jobs) job.input.clear();
    emit batchFinished(!m_failed && !m_stopping);
    if(m_closeWhenFinished) QDialog::reject();
}
bool ParallelDialog::writeSummary() {
    // Keep the live summary accurate even when the results file can no longer be written.
    QJsonArray jobs; for(const auto& job : m_jobs) jobs.append(job.record);
    m_summary["jobs"] = jobs;
    if(m_output.isEmpty() || m_persistFailed) return !m_persistFailed;
    QString error;
    if(!writeBytes(m_output + "/results.json", QJsonDocument(m_summary).toJson(), error)) {
        m_persistFailed = true; m_lastError = error; appendLog(error); return false;
    }
    return true;
}
void ParallelDialog::openSelectedResults() {
    const int row = m_table->currentRow();
    if(row >= 0 && row < m_jobs.size() && !m_jobs[row].folder.isEmpty()) {
        const QString report = m_jobs[row].folder + "/input.html";
        QDesktopServices::openUrl(QUrl::fromLocalFile(QFileInfo::exists(report) ? report : m_jobs[row].folder));
    } else if(!m_output.isEmpty()) QDesktopServices::openUrl(QUrl::fromLocalFile(m_output));
}
void ParallelDialog::reject() {
    if(m_running) { m_closeWhenFinished = true; cancelBatch(); }
    else QDialog::reject();
}
void ParallelDialog::closeEvent(QCloseEvent* event) {
    if(m_running) { event->ignore(); m_closeWhenFinished = true; cancelBatch(); }
    else QDialog::closeEvent(event);
}
