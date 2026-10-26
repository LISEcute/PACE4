#include <QApplication>
#include <QDir>
#include <QFile>
#include <QMetaObject>
#include <QTemporaryDir>
#include <QProgressDialog>
#include <QTimer>
#include <QFileInfo>
#include <QTextStream>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "tfusion.h"

QString FileArg;
int fontsizeGlobal = 9;
int useHighDpiScaling = 1;
bool PACESilentMode = true;
extern QString LISErootPATH, localPATH;
extern int _LMINN, _BatchMode, _BatchNsteps, _NCASC, _IDIST, _INPUT;
extern double _AJNUC;
extern FILE *f09, *f02, *f_particles;
extern QTextStream* htmlStream;
extern double _FYRST, _BARFACr, _ARATIO, _FACLA, _ALEVB[];

static void require(bool ok, const char* message)
{
    if(!ok) throw std::runtime_error(message);
}
static void near(double actual, double expected, const char* message)
{
    require(std::abs(actual - expected) < 1e-12, message);
}
static QStringList fields(const QString& line)
{
    return line.simplified().split(' ', Qt::SkipEmptyParts);
}

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    LISErootPATH = QCoreApplication::applicationDirPath();
    QTemporaryDir temp;
    require(temp.isValid(), "Temporary directory unavailable");
    localPATH = temp.path();
    try {
        for(int mode : {1, 2, 3, 5}) {
            FileArg = temp.filePath(QString("mode%1.pace").arg(mode));
            QFile fixture(FileArg);
            require(fixture.open(QIODevice::WriteOnly | QIODevice::Text), "Create fixture");
            QTextStream out(&fixture);
            out << "64 0 " << mode << " 0.625 0.75 1.125 0 10.25 0 0 0\n";
            if(mode == 1) {
                out << "20 48 50 124 0 0 0\n";
                out << "123.125 0 0 2.25 0.125 2 1 4.25\n";
                out << "\n"; // entrance optical potential
            } else {
                out << "28 60 30.125 0.25 " << (mode == 3 ? 3 : 10) << " 2\n";
                if(mode == 3) out << "3 0.2 0.3 0.5\n";
            }
            out << "\n\n\n";
            out.flush();
            fixture.close();
            {
                tfusion window;
                near(_ARATIO, 1.125, "Load level-density ratio");
                require(_LMINN == 2, "LMINN must load from its own field");
                require(QMetaObject::invokeMethod(&window, "on_actionSave_triggered", Qt::DirectConnection), "Invoke real Save action");
            }
            require(fixture.open(QIODevice::ReadOnly | QIODevice::Text), "Read saved fixture");
            QTextStream in(&fixture);
            auto header = fields(in.readLine());
            require(header.size() == 11, "Saved header field count");
            near(header[3].toDouble(), 0.625, "Save FYRST precision");
            near(header[4].toDouble(), 0.75, "Save barrier factor precision");
            near(header[5].toDouble(), 1.125, "Save level-density ratio precision");
            near(header[7].toDouble(), 10.25, "Save FACLA precision");
            auto nucleus = fields(in.readLine());
            if(mode == 1) {
                auto beam = fields(in.readLine());
                near(beam[0].toDouble(), 123.125, "Save loaded beam energy");
                near(beam[3].toDouble(), 2.25, "Save AGRAZ precision");
                near(beam[4].toDouble(), 0.125, "Save ELOSS precision");
                require(beam[6].toInt() == 1, "Save loaded quantum barrier mode");
                near(beam[7].toDouble(), 4.25, "Save loaded barrier curvature");
            } else {
                near(nucleus[2].toDouble(), 30.125, "Save loaded excitation energy");
                near(nucleus[3].toDouble(), 0.25, "Save recoil energy");
                require(nucleus[5].toInt() == 2, "Save LMINN");
            }
            fixture.close();
            {
                tfusion reloaded;
                near(_FYRST, 0.625, "Reload FYRST");
                near(_BARFACr, 0.75, "Reload barrier factor");
                near(_ARATIO, 1.125, "Reload level-density ratio");
                near(_FACLA, 10.25, "Reload FACLA");
                require(_LMINN == 2, "Reload LMINN");
                if(mode == 3) {
                    near(_ALEVB[1], 0.2, "Reload spin weight 1");
                    near(_ALEVB[2], 0.3, "Reload spin weight 2");
                    near(_ALEVB[3], 0.5, "Reload spin weight 3");
                }
            }
            std::cout << "PASS mode " << mode << " load/save/reload\n";
        }
        {
            tfusion engine;
            _BatchMode = 1;
            _BatchNsteps = 1;
            require(QMetaObject::invokeMethod(&engine, "on_actionExecute_triggered", Qt::DirectConnection), "Invoke Execute with invalid batch");
            require(!QFileInfo::exists(engine.FFileNameOut), "Invalid batch must not create outputs");
            require(!engine.runSilent(), "Invalid batch must return failure");
            _BatchMode = 0;
            std::cout << "PASS one-step batch rejection through Execute and silent mode\n";

            auto outputsClosed = []() {
                require(!f09 && !f02 && !f_particles && !htmlStream, "All output resources must be released");
            };
            auto calculate = [&]() {
                return engine.p->PACE3(engine.FFileNameOut.toLocal8Bit().constData(),
                                       engine.FFileNameEve.toLocal8Bit().constData(),
                                       engine.FFileNameCS.toLocal8Bit().constData(),
                                       engine.FFileNameParticles.toLocal8Bit().constData(),
                                       engine.FFileNameHtml.toLocal8Bit().constData());
            };
            _INPUT = 5;
            _AJNUC = 200;
            require(calculate() < 0, "COMPOS failure must remain a failure");
            outputsClosed();
            std::cout << "PASS COMPOS error propagation and output cleanup\n";

            _AJNUC = 10;
            _NCASC = 10000;
            _IDIST = 1;
            bool cancellationRequested = false;
            QTimer::singleShot(0, [&]() {
                for(QWidget* widget : QApplication::topLevelWidgets()) {
                    if(auto progress = qobject_cast<QProgressDialog*>(widget)) {
                        progress->cancel();
                        cancellationRequested = true;
                    }
                }
            });
            PACESilentMode = false;
            const int canceled = calculate();
            PACESilentMode = true;
            require(cancellationRequested && canceled == 10, "Simulation cancellation must preserve canceled status");
            outputsClosed();
            require(QFile::remove(engine.FFileNameEve), "Canceled event file must be closed");
            std::cout << "PASS actual cancellation and output cleanup\n";
        }
    } catch(const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
