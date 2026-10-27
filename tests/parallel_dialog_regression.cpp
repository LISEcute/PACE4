// Native QProcess integration tests. No Python required.
// Optional real physics: set PACE_DIALOG_REAL_EXE to a deployed PACE4.exe.
#include "parallel_dialog.h"
#include <QApplication>
#include <QTest>
#include <QTemporaryDir>
#include <QFile>
#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QSpinBox>
#include <QPushButton>
#include <QTableWidget>
#include <QProgressBar>
#include <QThread>
#include <QElapsedTimer>
#include <QSet>
#include <cmath>

static QByteArray read(const QString& path) {
    QFile f(path); if (!f.open(QIODevice::ReadOnly)) return {}; return f.readAll();
}
static bool write(const QString& path, const QByteArray& bytes) {
    QFile f(path); return f.open(QIODevice::WriteOnly) && f.write(bytes)==bytes.size();
}
static QByteArray fixture(int idist=1, int events=128, const QByteArray& control={}) {
    return QByteArray::number(events)+" "+QByteArray::number(idist)+" 5 1 -1 1 0 10 0 1 31\n28 60 30.125 0 20 0\n\n\n\n"+control;
}
// This very executable doubles as a controlled independent child.
static int worker(const QString& path) {
    const QByteArray data=read(path);
    if(data.contains("FAIL_WORKER")) return 7;
    if(data.contains("SLEEP_WORKER")) QThread::sleep(30);
    else QThread::msleep(180);
    if(data.contains("MISSING_OUTPUT")) return 0;
    if(qEnvironmentVariableIsSet("PACE_BENCHMARK_SPINS")) return 9;
    const QByteArray payload=data+"\nseed="+qgetenv("PACE4_SEED");
    const QString base=QFileInfo(path).absolutePath()+"/"+QFileInfo(path).completeBaseName();
    for(const auto& suffix:{".cs4",".html",".rtf"}) if(!write(base+suffix,payload)) return 10;
    QDir().mkpath("runtime/files");
    return write("runtime/files/PACE4evt.dat",payload)?0:11;
}

class DialogTests: public QObject {
    Q_OBJECT
    QTemporaryDir temp;
    int sequence=0;
    QString source(const QByteArray& bytes) {
        const QString p=temp.filePath(QString("source-%1.pace").arg(++sequence));
        if(!write(p,bytes)) qFatal("Cannot create fixture");
        return p;
    }
    QString output() { return temp.filePath(QString("result-%1").arg(++sequence)); }
    void populate(ParallelDialog& d, const QStringList& files, int workers) {
        d.addInputFiles(files); d.setOutputDirectory(output()); d.setWorkerCount(workers);
        for(int i=0;i<files.size();++i) d.setJobSeed(i,104729+i*130363);
    }
    bool wait(ParallelDialog& d) {
        QElapsedTimer timer; timer.start();
        while(d.isRunning() && timer.elapsed()<90000) QTest::qWait(20);
        return !d.isRunning();
    }
    static QJsonArray jobs(const ParallelDialog& d) {return d.summary().value("jobs").toArray();}
    static QString folder(const QJsonValue& j) {return j.toObject().value("folder").toString();}
private slots:
    void initTestCase() {
        QVERIFY(temp.isValid());
        const QString cfg=QCoreApplication::applicationDirPath()+"/lisecfg";
        QVERIFY(QDir().mkpath(cfg));
        QVERIFY(write(cfg+"/AME_DB.sqlite","fake independent database"));
    }
    void controlsAndScreenshot() {
        ParallelDialog d;
        populate(d,{source(fixture()),source(fixture(2))},2);
        auto table=d.findChild<QTableWidget*>("inputTable"); QVERIFY(table); QCOMPARE(table->rowCount(),2);
        auto workers=d.findChild<QSpinBox*>("workerCount"); QVERIFY(workers); QVERIFY(workers->minimum()>=1);
        d.setWorkerCount(0); QVERIFY(workers->value()>=1);
        d.setWorkerCount(2);
        for(const auto& name:{"addInputs","duplicateInputs","removeInputs","startBatch","cancelBatch","closeDialog"})
            QVERIFY2(d.findChild<QPushButton*>(name),name);
        QVERIFY(d.findChild<QProgressBar*>("progressBar"));
        d.show(); QTest::qWait(150);
        QVERIFY(d.grab().save(QCoreApplication::applicationDirPath()+"/parallel-dialog.png"));
        table->selectRow(0);
        d.findChild<QPushButton*>("duplicateInputs")->click(); QCOMPARE(table->rowCount(),3);
        table->selectRow(2);
        d.findChild<QPushButton*>("removeInputs")->click(); QCOMPARE(table->rowCount(),2);
    }
    void validationAndNoOverwrite() {
        ParallelDialog empty; empty.setOutputDirectory(output()); QVERIFY(!empty.startBatch()); QVERIFY(!empty.lastError().isEmpty());
        const QString input=source(fixture());
        ParallelDialog duplicate; populate(duplicate,{input,input},2);
        duplicate.setJobSeed(0,123); duplicate.setJobSeed(1,123);
        QVERIFY(!duplicate.startBatch()); QVERIFY(!duplicate.isRunning());
        duplicate.setJobSeed(0,123459876); duplicate.setJobSeed(1,123459877);
        QVERIFY(!duplicate.startBatch()); QVERIFY(!duplicate.isRunning());
        const QString existing=output(); QVERIFY(QDir().mkpath(existing)); QVERIFY(write(existing+"/keep.txt","keep"));
        ParallelDialog overwrite; populate(overwrite,{input},1); overwrite.setOutputDirectory(existing);
        QVERIFY(!overwrite.startBatch()); QCOMPARE(read(existing+"/keep.txt"),QByteArray("keep"));
    }
    void serialParallelPreservesInputsAndCutoff() {
        const QStringList inputs={source(fixture(0,100000)),source(fixture(1,100000)),source(fixture(2,100000))};
        ParallelDialog serial, parallel; populate(serial,inputs,1); populate(parallel,inputs,2);
        QSet<int> active; int peak=0;
        connect(&parallel,&ParallelDialog::jobStateChanged,this,[&](int row,const QString& state){
            if(state=="running") {active.insert(row); peak=qMax(peak,int(active.size()));}
            else active.remove(row);
        });
        qputenv("PACE_BENCHMARK_SPINS","must-be-removed");
        QVERIFY2(serial.startBatch(),qPrintable(serial.lastError())); QVERIFY(wait(serial));
        QVERIFY2(parallel.startBatch(),qPrintable(parallel.lastError())); QVERIFY(wait(parallel));
        qunsetenv("PACE_BENCHMARK_SPINS");
        QCOMPARE(serial.summary().value("status").toString(),QString("succeeded"));
        QCOMPARE(parallel.summary().value("status").toString(),QString("succeeded"));
        QCOMPARE(peak,2); QVERIFY(active.isEmpty());
        auto a=jobs(serial),b=jobs(parallel); QCOMPARE(a.size(),3); QCOMPARE(b.size(),3);
        QSet<QString> folders;
        for(int i=0;i<3;++i) {
            const auto j=b[i].toObject(); folders.insert(folder(b[i]));
            QCOMPARE(read(folder(b[i])+"/input.pace"),read(inputs[i]));
            QCOMPARE(j.value("ncasc").toInt(),100000);
            QVERIFY(std::abs(j.value("effective_discriminator").toDouble()-1/(3*std::sqrt(100000.)+10))<1e-12);
            QCOMPARE(read(folder(a[i])+"/input.cs4"),read(folder(b[i])+"/input.cs4"));
            QCOMPARE(read(folder(a[i])+"/runtime/files/PACE4evt.dat"),read(folder(b[i])+"/runtime/files/PACE4evt.dat"));
            QVERIFY(!read(folder(b[i])+"/input.cs4").isEmpty());
        }
        QCOMPARE(folders.size(),3);
    }
    void failureStopsQueue() {
        ParallelDialog d; populate(d,{source(fixture(1,128,"FAIL_WORKER")),source(fixture(1,128,"SLEEP_WORKER")),source(fixture())},2);
        QVERIFY(d.startBatch()); QVERIFY(wait(d));
        QCOMPARE(jobs(d).size(),3);
        QCOMPARE(jobs(d)[0].toObject().value("status").toString(),QString("failed"));
        QCOMPARE(jobs(d)[1].toObject().value("status").toString(),QString("canceled"));
        QCOMPARE(jobs(d)[2].toObject().value("status").toString(),QString("canceled"));
        QVERIFY(d.findChild<QPushButton*>("closeDialog")->isEnabled());
    }
    void missingOutputsAndLaunchError() {
        ParallelDialog missing; populate(missing,{source(fixture(1,128,"MISSING_OUTPUT"))},1);
        QVERIFY(missing.startBatch()); QVERIFY(wait(missing));
        QCOMPARE(jobs(missing)[0].toObject().value("status").toString(),QString("failed"));
        ParallelDialog bad(nullptr,temp.filePath("not-an-executable.exe")); populate(bad,{source(fixture())},1);
        QVERIFY(!bad.startBatch()); QVERIFY(!bad.isRunning()); QVERIFY(!bad.lastError().isEmpty());
    }
    void cancelReapsChildren() {
        ParallelDialog d; populate(d,{source(fixture(1,128,"SLEEP_WORKER")),source(fixture(1,128,"SLEEP_WORKER")),source(fixture())},2);
        QVERIFY(d.startBatch()); QTest::qWait(250); QVERIFY(d.isRunning());
        d.cancelBatch(); QVERIFY(wait(d));
        for(const auto& j:jobs(d)) QCOMPARE(j.toObject().value("status").toString(),QString("canceled"));
        QVERIFY(d.findChild<QPushButton*>("closeDialog")->isEnabled());
        // Children are reaped; runtime data can be removed immediately on Windows.
        for(const auto& j:jobs(d)) {
            QVERIFY(QDir::cleanPath(folder(j)).startsWith(QDir::cleanPath(temp.path())+"/"));
            QVERIFY(QDir(folder(j)).removeRecursively());
        }
    }
    void persistenceFailureKeepsTerminalStates() {
        ParallelDialog d; populate(d,{source(fixture(1,128,"SLEEP_WORKER")),source(fixture(1,128,"SLEEP_WORKER"))},2);
        QVERIFY(d.startBatch()); QTest::qWait(250);
        const QString result=QFileInfo(folder(jobs(d)[0])).dir().absolutePath()+"/../results.json";
        QVERIFY(QDir::cleanPath(result).startsWith(QDir::cleanPath(temp.path())+"/"));
        QVERIFY(QFile::remove(result)); QVERIFY(QDir().mkpath(result));
        d.cancelBatch(); QVERIFY(wait(d));
        QCOMPARE(d.summary().value("status").toString(),QString("failed"));
        for(const auto& j:jobs(d)) QVERIFY(j.toObject().value("status").toString()!="running");
    }
    void timeoutStopsQueue() {
        ParallelDialog d; populate(d,{source(fixture(1,128,"SLEEP_WORKER")),source(fixture())},1);
        d.setTimeoutSeconds(1); QVERIFY(d.startBatch()); QVERIFY(wait(d));
        QCOMPARE(jobs(d)[0].toObject().value("status").toString(),QString("timed_out"));
        QCOMPARE(jobs(d)[1].toObject().value("status").toString(),QString("canceled"));
    }
    void closeStopsChildren() {
        ParallelDialog d; populate(d,{source(fixture(1,128,"SLEEP_WORKER"))},1);
        d.show(); QVERIFY(d.startBatch()); QTest::qWait(100); d.close(); QVERIFY(wait(d));
        QCOMPARE(jobs(d)[0].toObject().value("status").toString(),QString("canceled"));
    }
    void realPhysicsSerialParallel() {
        const QString exe=qEnvironmentVariable("PACE_DIALOG_REAL_EXE");
        if(exe.isEmpty()) QSKIP("Set PACE_DIALOG_REAL_EXE for deployed physics equivalence");
        QVERIFY(QFileInfo::exists(exe));
        const QStringList inputs={source(fixture(0,40000)),source(fixture(1)),source(fixture(2))};
        ParallelDialog serial(nullptr,exe),parallel(nullptr,exe); populate(serial,inputs,1); populate(parallel,inputs,2);
        QVERIFY2(serial.startBatch(),qPrintable(serial.lastError())); QVERIFY(wait(serial));
        QVERIFY2(parallel.startBatch(),qPrintable(parallel.lastError())); QVERIFY(wait(parallel));
        QVERIFY(write(QCoreApplication::applicationDirPath()+"/real-parallel-summary.json",QJsonDocument(parallel.summary()).toJson()));
        QVERIFY2(parallel.summary().value("status").toString()=="succeeded",qPrintable(parallel.lastError()));
        QVERIFY(std::abs(jobs(parallel)[0].toObject().value("effective_discriminator").toDouble()-1./610)<1e-12);
        QCOMPARE(read(folder(jobs(parallel)[0])+"/input.pace"),read(inputs[0]));
        for(int i=0;i<3;++i) {
            QCOMPARE(jobs(serial)[i].toObject().value("status").toString(),QString("succeeded"));
            QCOMPARE(jobs(parallel)[i].toObject().value("status").toString(),QString("succeeded"));
            for(const auto& name:{"/input.cs4","/input.particles","/runtime/files/PACE4evt.dat"})
                QCOMPARE(read(folder(jobs(serial)[i])+name),read(folder(jobs(parallel)[i])+name));
            QVERIFY(!read(folder(jobs(parallel)[i])+"/input.cs4").isEmpty());
            if(i>0) {
                QVERIFY(!read(folder(jobs(parallel)[i])+"/runtime/files/PACE4evt.dat").isEmpty());
                QVERIFY(!read(folder(jobs(parallel)[i])+"/input.particles").isEmpty());
            }
        }
    }
};
int main(int argc,char** argv) {
    if(argc>1 && QString::fromLocal8Bit(argv[1]).endsWith(".pace")) {
        QCoreApplication app(argc,argv); return worker(QString::fromLocal8Bit(argv[1]));
    }
    QApplication app(argc,argv); DialogTests tests; return QTest::qExec(&tests,argc,argv);
}
#include "parallel_dialog_regression.moc"






