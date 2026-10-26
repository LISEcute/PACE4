#include "tfusion.h"

#include <QApplication>
//#include <QDesktopWidget>
#include <QStyleFactory>
#include <QDir>
#include <QSettings>
#include <QStandardPaths>
#include <QCoreApplication>


extern void getInitialDir(void);
extern void getDPIscaling(void);

extern QString LISErootPATH;
extern QString MyDocCompPATH;
extern const char *LISEini;

QString FileArg="";
int fontsizeGlobal=9;
int useHighDpiScaling=1;
bool PACESilentMode=false;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int main(int argc, char *argv[])
{

//getDPIscaling();
//if(!useHighDpiScaling)
//  QCoreApplication::setAttribute(Qt::AA_Use96Dpi);

//   QApplication::setAttribute(useHighDpiScaling ? Qt::AA_EnableHighDpiScaling :
//                                                   Qt::AA_DisableHighDpiScaling);
//--------------------------------------------------------------------------------------
    QApplication app(argc, argv);

    for(int i=1; i<argc; ++i)
      {
        const QString arg = QDir::fromNativeSeparators(QString::fromLocal8Bit(argv[i])).trimmed();
        if(arg.compare(QStringLiteral("-s"), Qt::CaseInsensitive) == 0 ||
           arg.compare(QStringLiteral("/s"), Qt::CaseInsensitive) == 0)
          {
            PACESilentMode = true;
          }
        else if(FileArg.isEmpty())
          {
            FileArg = arg;
          }
      }
    LISErootPATH = QCoreApplication::applicationDirPath();
    if(LISErootPATH.size()<=0)
        {
        QString arg0 = QDir::fromNativeSeparators(argv[0]);
        QDir dir(arg0);
        LISErootPATH = dir.currentPath();
        }

    //---------------------------------------------------------------------------------- paths + read font size start
    getInitialDir();
    //---------------------------------------------------------------------------------- paths + read font size stop

    //------------------------------------------------------------------ style start

    app.setFont(QFont("Arial", fontsizeGlobal, QFont::Normal));
    app.setOrganizationName("FRIB/MSU");
    app.setStyle(QStyleFactory::create("Fusion"));
    QFile FileStyle(":/styles/mainstyle.qss");
    if(FileStyle.open(QFile::ReadOnly)) {
        qApp->setStyleSheet(QString::fromUtf8(FileStyle.readAll()));
    } else {
        qWarning() << "PACE4: could not open bundled stylesheet:" << FileStyle.fileName();
    }

    int fontsize = 9;

    app.setFont(QFont("Arial", fontsize, QFont::Normal));
    //------------------------------------------------------------------ style end

    //---------------------------------------------------------------- Daniel's begin; OT modified

    QPalette palette;
    palette.setColor(QPalette::Window, QColor(240, 245, 240)); // set window background color
    palette.setColor(QPalette::WindowText, QColor(0, 0, 0)); // set text color
    palette.setColor(QPalette::Button, QColor(235, 240, 235)); // set button background color
    palette.setColor(QPalette::ButtonText, QColor(0, 10, 0)); // set button text color
    palette.setColor(QPalette::Text, QColor(0, 0, 0)); // set text color for QLineEdit
    palette.setColor(QPalette::Base, QColor(254, 255, 254));
    palette.setColor(QPalette::Light, QColor(255, 255, 255)); //set the light part of the sunken border palette
    palette.setColor(QPalette::Dark, QColor(159, 159, 159)); //set the dark part of the sunken border palette

    // Set disabled line edit style
    /*  QString existingStyleSheet = qtApp.styleSheet();
  QString disabledLineEditStyle = QString("QLineEdit:disabled { color: %1; }").arg(LE_tc_disabled.name());
  QString disabledButtonStyle = QString("QPushButton:disabled { background-color: %1; color: %2; }")
      .arg(LE_bg.name())
      .arg(LE_tc_disabled.name());
  QString updatedStyleSheet = existingStyleSheet + disabledLineEditStyle + disabledButtonStyle;
  qtApp.setStyleSheet(updatedStyleSheet);
*/

    QColor bg_disabled(215, 220, 215);
    QColor tc_disabled(145,150,145);

    palette.setColor(QPalette::Disabled, QPalette::Button, bg_disabled); // Set disabled QLineEdit background color
    //  palette.setColor(QPalette::Disabled, QPalette::Base, bg_disabled); // Set disabled QLineEdit background color
    palette.setColor(QPalette::Disabled, QPalette::Text, tc_disabled); // Set disabled QLineEdit text color
    palette.setColor(QPalette::Disabled, QPalette::WindowText, tc_disabled); // Set disabled QLineEdit text color
    palette.setColor(QPalette::Disabled,  QPalette::ButtonText, tc_disabled); // set button text color
    app.setPalette(palette);

     //------------------------------------------------------------------ style end

    if(PACESilentMode && FileArg.isEmpty())
      {
        qWarning() << "PACE4 silent mode requires an input .pace file.";
        return 2;
      }

    tfusion widget;
    widget.setMinimumWidth(750);
    widget.setMinimumHeight(754);

    if(PACESilentMode)
      return widget.runSilent() ? 0 : 1;

    widget.show();

    return app.exec();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
