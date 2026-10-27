#ifndef PARALLEL_DIALOG_H
#define PARALLEL_DIALOG_H

#include <QDialog>
#include <QElapsedTimer>
#include <QJsonObject>
#include <QList>
#include <QPointer>
#include <QStringList>

class QCloseEvent;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QProcess;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTableWidget;
class QTimer;

class ParallelDialog : public QDialog
{
    Q_OBJECT
public:
    explicit ParallelDialog(QWidget* parent = nullptr, const QString& executable = QString());
    ~ParallelDialog() override;
    void addInputFiles(const QStringList& paths);
    void setOutputDirectory(const QString& path);
    void setWorkerCount(int count);
    void setJobSeed(int row, quint32 seed);
    void setTimeoutSeconds(int seconds);
    bool startBatch();
    void cancelBatch();
    bool isRunning() const { return m_running; }
    QString lastError() const { return m_lastError; }
    QJsonObject summary() const { return m_summary; }

signals:
    void batchFinished(bool success);
    void jobStateChanged(int row, const QString& state);
    void validationFailed(const QString& message);

protected:
    void closeEvent(QCloseEvent* event) override;
    void reject() override;

private:
    struct Job {
        QString id, source, folder, state = QStringLiteral("pending");
        quint32 seed = 1;
        int ncasc = 0;
        QByteArray input;
        QJsonObject record;
        QPointer<QProcess> process;
        QElapsedTimer elapsed;
    };
    void refreshRows();
    void setState(int row, const QString& state);
    void setEditingEnabled(bool enabled);
    void showError(const QString& message);
    void appendLog(const QString& message);
    void schedule();
    void launch(int row);
    void processFinished(int row, int exitCode, bool crashed);
    void stopBatch(const QString& reason, bool failure);
    void completeIfStopped();
    void finishBatch();
    bool writeSummary();
    void collectOutputs(Job& job);
    void openSelectedResults();
    void checkTimeouts();
    QString proposedOutput(const QString& parent) const;

    QList<Job> m_jobs;
    QString m_executable, m_output, m_lastError;
    QByteArray m_binary, m_database;
    QJsonObject m_summary;
    QElapsedTimer m_totalClock, m_executionClock;
    quint64 m_nextIdentity = 1;
    quint32 m_nextSeed = 436;
    int m_nextJob = 0, m_active = 0;
    bool m_running = false, m_stopping = false, m_failed = false;
    bool m_closeWhenFinished = false, m_persistFailed = false;
    QTableWidget* m_table;
    QSpinBox* m_workers;
    QSpinBox* m_timeout;
    QLineEdit* m_outputEdit;
    QLabel* m_status;
    QPlainTextEdit* m_log;
    QProgressBar* m_progress;
    QPushButton *m_add, *m_duplicate, *m_remove, *m_browse, *m_start, *m_cancel, *m_open, *m_close;
    QTimer* m_timer;
};
#endif
