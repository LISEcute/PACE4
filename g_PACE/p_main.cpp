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
    LISErootPATH = QCoreApplication::applicationDirPath();
    if(LISErootPATH.size()<=0)
        {
        QString arg0 = QDir::fromNativeSeparators(argv[0]);
        QDir dir(arg0);
        LISErootPATH = dir.currentPath();
        }

    if(argc>1) FileArg = QDir::fromNativeSeparators(argv[1]);

    //---------------------------------------------------------------------------------- paths + read font size start
    getInitialDir();
    //---------------------------------------------------------------------------------- paths + read font size stop

    //------------------------------------------------------------------ style start

    app.setFont(QFont("Arial", fontsizeGlobal, QFont::Normal));
    app.setOrganizationName("FRIB/MSU");
    app.setStyle(QStyleFactory::create("Fusion"));
    QFile FileStyle(":/w_Main/mainstyle.qss");
    FileStyle.open(QFile::ReadOnly);
    QString StyleSheet = QLatin1String(FileStyle.readAll());
    qApp->setStyleSheet(StyleSheet);

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

    tfusion widget;
    widget.setMinimumWidth(750);
    widget.setMinimumHeight(754);
    widget.show();

    return app.exec();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
