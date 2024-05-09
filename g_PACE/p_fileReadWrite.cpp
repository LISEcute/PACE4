#include "tfusion.h"
#include <QDir>
#include <QStandardPaths>
#include <iostream>
#include <QSettings>
#include <QDebug>
#include <QFileInfo>

const char *dir_files="/files";
const char *FileNameAbsent = "/pUntitled";
const char *LISEini="/lisepp.ini";
QString LISErootPATH;
QString MyDocCompPATH;
QString localPATH;
QString basePATH;

//QString globalPATH;

extern int fontsizeGlobal;
extern int useHighDpiScaling;
FILE *mfopen(const QString& filename, const char* operand);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getDPIscaling(void)
{
      MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();
      QString FN1 = MyDocCompPATH + "/LISEcute";  FN1 += LISEini;
      QSettings myLiseIni1(FN1,QSettings::IniFormat);
      useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void getInitialDir(void)
{
    //---------------------------------------------------------  paths begin

        MyDocCompPATH = QStandardPaths::standardLocations(QStandardPaths::DocumentsLocation).constFirst();
        //qDebug() << MyDocCompPATH;
        QString FileCheck(LISErootPATH); FileCheck += LISEini;
        FILE *fcheck=mfopen(FileCheck,"at");
        int work_in_LISEroot_main = 0;

        if(fcheck) {                    // work in root directory
              fclose(fcheck);
              QSettings myLiseIni0(FileCheck,QSettings::IniFormat);
              work_in_LISEroot_main = myLiseIni0.value("Version/WorkInROOT",0).toInt();
              if(work_in_LISEroot_main) localPATH = LISErootPATH;
              }

        if(work_in_LISEroot_main==0)
              {
              localPATH = MyDocCompPATH;
              localPATH += "/LISEcute";
              }

     basePATH = localPATH + "/lisecfg/";

     //--------------------------------------------------------- lise.ini  begin
       QString FN1=localPATH;  FN1 += LISEini;
       QSettings myLiseIni1(FN1,QSettings::IniFormat);
       fontsizeGlobal     = myLiseIni1.value("font/size",     fontsizeGlobal).toInt();
       useHighDpiScaling  = myLiseIni1.value("font/scaling",  useHighDpiScaling).toInt();
      //--------------------------------------------------------- lise.ini  end

     localPATH += dir_files;
     QDir pathDir(localPATH);
     if(!pathDir.exists()) pathDir.mkdir(localPATH);

    //---------------------------------------------------------  paths end

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void tfusion::setFileName(const QString& fileName)
{    
    //    QString windowName = FFileName.split("/",Qt::SkipEmptyParts).last();

    FFileName = fileName;
    QFileInfo fI(FFileName);
    QString windowName = fI.baseName();
    if(!windowName.contains(&FileNameAbsent[1]))  this->setWindowTitle("PACE4 - " + windowName);
    else                                          this->setWindowTitle("PACE4");

    // to avoid problems with //intranet.fff.fff.ff/ path
    int posDot   =  FFileName.lastIndexOf('.');
    int posSlash =  FFileName.lastIndexOf('/');
    QString base;

//    if(posDot > 0 && posSlash < posDot) base = FFileName.split(".",Qt::SkipEmptyParts).at(0);
    if(posDot > 0 && posSlash < posDot) base = FFileName.left(posDot);
    else                                base = FFileName;

//    qDebug() << "FFileName & base " << FFileName << base << windowName;
    FFileNameOut        = base + ".rtf";
    FFileNameParticles  = base + ".particles";
    FFileNameCS         = base + ".cs4";
    FFileNameHtml       = base + ".html";
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*FILE *mfopen(const QString& filename, const char* operand)
{

#if !defined(__CYGWIN__) && !defined(_WIN32) && !defined(_WIN64)

  FILE *f =   fopen(filename.toStdString().c_str(),operand);

#else

  QString woperand(operand);
  FILE *f = _wfopen(filename.toStdWString().c_str(), woperand.toStdWString().c_str() );

#endif

return f;
}
*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
