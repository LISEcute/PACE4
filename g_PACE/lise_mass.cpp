#include <QFile>
#include <QDebug>
#include <iostream>
#include <QFileInfo>
#include <QMessageBox>
#include <QApplication>
#include <stdlib.h>
#include <QString>
#include <QtMath>

#include <QDebug>
#include <QDir>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlRecord>
#include <QVector>

#include "g_PACE/Constant.h"
#include "tfusion.h"


#define pow(x,y)  exp((y)*log(x))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))


double W_ME(int A, int Z);

double  receiveQ(int Z, int A, int &opt);   // opt=0   default from database
bool FinderPace(int A, int Z, std::vector<std::vector<QVariant>>& Array );

//void close_dbf( );
bool init_base(QString name);
double Weitseker(int A, int Z);
double MASSES(int IZ, int IN, int &opt);   // opt - gde bylo prochitano
double MASSES1(int IZ, int IN);
double  YMASS(int Z, int A);

//bool found=false;
//DBF *s_DBF=nullptr;
extern QString basePATH;
extern QString basePATH2;
extern QString LISErootPATH;

extern QSqlDatabase db;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
bool tfusion::init_base()
{
  QString AME_DbName = "AME_DB.sqlite"; // Name of the AME Database File
  QString AMEDbPath = LISErootPATH + "/lisecfg" + "/" + AME_DbName; // File Path of the AME Database file

//  qDebug() << AMEDbPath;
//#define ACCESS "DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};"
//  "FIL={MS Access};DBQ=" + AMEDbPath

  db= QSqlDatabase::addDatabase("QSQLITE", "AME2016_connection");
  db.setDatabaseName(AMEDbPath);

//  db = QSqlDatabase::addDatabase("QODBC");
//  db.setDatabaseName(ACCESS);
  db.open();

  if(db.isOpenError())
    QMessageBox::warning(this, "Database file Didn't open",AMEDbPath+" Can't be located");

  return db.isOpen();
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//========================== Search Function =============================
//========================================================================
bool FinderPace(int A, int Z, std::vector<std::vector<QVariant>>& Array )
{
  QSqlQuery query(db);
  query.prepare("SELECT * FROM AME2016 WHERE A = ? AND Z = ? ORDER BY [INDEX] ASC");
  query.addBindValue(A);
  query.addBindValue(Z);

  //====================================================================================
  if (!query.exec()) { qDebug() << "Search Didn't work" << A << Z;
                     return false;  }
  //====================================================================================
  while (query.next())
    {
      std::vector<QVariant> row; // this creates a 1d array

      for (int i= 0; i < query.record().count();i++)
        row.push_back(query.value(i)); // here we add the values into the 1d vecotr

      Array.push_back(row); // we append the 1d vector into the 2d vector
    }

  return Array.size() > 0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double  receiveQ(int Z, int A, int &opt)   // opt=0   default from database
{                                         // opt=1   calculation
  if(Z==0) return A*ME_n;

  std::vector<std::vector<QVariant>> Array;
  bool found = FinderPace(A,Z,Array);

  if(opt==0 && found)
    if(Array[0][5].toString()!='*')
      return Array[0][5].toDouble();

  opt=1;

  return W_ME(A,Z);
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double W_ME(int A, int Z)
{
  double dW0 = Weitseker(A,Z);
  double result = Z*ME_p + (A-Z)*ME_n - dW0;
  return result;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double Weitseker(int A, int Z)
{
  if(A*Z <= 0) return -1000;

  double alpha=15.5; //15.75;
  double beta =16.8; //17.8;
  double gamma=0.72; //0.71;
  double epsilon=23; //94.8/4;
  double delta=34.0;

  double A1=alpha*A;
  double A2=-beta*pow(A,2./3.);
  double A3=-gamma*Z*Z*pow(A,-1./3.);
  double t= A - 2.*Z;
  double A4=-epsilon*t*t/A;

  double d=0;
  if( ((int(A)+1)%2)  *
      ((int(Z)+1)%2)  == 1) d=delta;
  if( ((int(A)  )%2)  *
      ((int(Z)  )%2)  == 1) d=-delta;
  double A5=0;
  if(d!=0) A5=d*pow(A,-0.75);
  double result=A1+A2+A3+A4+A5;
  return result;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double MASSES(int IZ, int IN, int &opt)
{
  extern int _NOSHL;
  opt = _NOSHL;

  return (_NOSHL==0 ? receiveQ(IZ, IZ+IN, opt) : MASSES1(IZ, IN));
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double MASSES1(int IZ, int IN)
{
  extern FILE *f09;

  double DEFEC=YMASS(IZ,IZ+IN);

  if(DEFEC==999.)fprintf(f09,
                         "\\par  Warning - no physical mass defect available. NOSHL = 1. try to use NOSHL=0"
                         "\\par  Decay mode to this nucleus  and following chain suppressed"
                         "\\par  Edit mass table to correct for missing mass");

  return DEFEC;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

double  YMASS(int Z, int A)
{
  //     FROM PROGRAM ALICE (F. PLASIL & M. BLANN)
  // common/XQSHL/NOSHL;
  extern int _NOSHL;
  //extern FILE *f09;

  double EM[11],XK[11],Y[3],F[3],EMP[11];

  //double ZEE=Z;
  double YMASS=999.;

  if(Z>=21.) {

      EM[1]=0.0;
      EM[2]=2.0;
      EM[3]=8.0;
      EM[4]=14.0;
      EM[5]=28.0;
      EM[6]=50.0;
      EM[7]=82.0;
      EM[8]=126.0;
      EM[9]=184.0;
      EM[10]=258.0;
      double CAY1=1.15303;
      double CAY2=0.0;
      double CAY3=200.0;
      double CAY4=11.0;
      double CAY5=8.07144;
      double CAY6=7.28899;
      double GAMMA=1.7826;
      double A1=15.4941;
      double A2=17.9439;
      double A3=0.7053;
      double D=0.444;
      double C=5.8;
      double SMALC=0.325;
      int I, J, IPQ;

      for(I=1; I<=10; I++)  EMP[I]=(EM[I]>0 ? pow(EM[I],(5./3.)) : 0);
      for(I=1; I<=9; I++)   XK[I]=0.6*(EMP[I+1]-EMP[I])/(EM[I+1]-EM[I]);

      //double RZ=0.863987/A3;
      //int L=0;
      double UN=A-Z;

      double A3RT=pow(A,(1./3.));
      //   double A2RT=sqrt(A);
      double A3RT2=pow(A3RT,2.0);
      double ZSQ=Z*Z;

      double SYM=pow(((UN-Z)/A),2);
      double ACOR=1.0-GAMMA*SYM;
      double PARMAS=CAY5*UN+CAY6*Z;
      double VOLNUC=-1.0*A1*ACOR*A;
      double SUFNUC=A2*ACOR*A3RT2;
      double COULMB=A3*ZSQ/A3RT;
      double FUZSUR=-1.0*CAY1*ZSQ/A;
      int    N=UN+.01;
      int    IZ=Z+.01;
      double ODDEV=-1.0*(1.0+2.0*(N/2)-UN+2.*(IZ/2)-Z)/sqrt(A)*CAY4;
      double WTERM;

      if(SYM>0.4)WTERM=0.;
      else       WTERM=-1.*CAY2*A3RT2*exp(-1.*CAY3*SYM);

      double WOTNUC=PARMAS+COULMB+FUZSUR+ODDEV+WTERM;
      double SMASS=WOTNUC+VOLNUC+SUFNUC;
      double C2=(SUFNUC+WTERM)/(pow(A,(2./3.)));
      double X=COULMB/(2.0*(SUFNUC+WTERM));
      //   double BARR=0.0;
      Y[1]=UN;
      Y[2]=Z;


      for(J=1; J<=2; J++) {
          for(I=1; I<=9; I++) if (Y[J]-EM[I+1] <=0) break;

          if(I>9) goto fin;
          F[J]=XK[I]*(Y[J]-EM[I])-0.6*(pow(Y[J],(5./3.))-EMP[I]);
        }


      double S=pow((2.0/A),(2.0/3.0))*(F[1]+F[2])-SMALC*pow(A,(1./3.));
      double EE=2.*C2*pow(D,2)*(1.0-X);
      double FF=0.42591771*C2*pow(D,3)*(1.+2.*X)/A3RT;
      double SSHELL=C*S;
      double V=SSHELL/EE;
      double EPS=1.5*FF/EE;
      //double QCALC=0.0,
      // double THETA=0.0,
      double SHLL=0, TO=0, T=0, GL=0;
      //   double ALPHA0;//, ALPHA;//, SIGMA;

      //-----------------
      if( EE*(1.-3.*V) > 0) SHLL=SSHELL;
      //-----------------
      else
        {

          TO=1.0;

L175:

          for(IPQ=1; IPQ<=10; IPQ++)
            {

              if (TO<=12.)  T=TO-(1.-EPS*TO-V*(3.-2.*TO*TO)*exp(-TO*TO))/
                  (-EPS+V*(10.*TO-4. *pow(TO,3))*exp(-TO*TO));

              else          T=TO-(1.-EPS*TO)/(-EPS);

              if (T<=0.0)               goto L190;
              if (fabs(T-TO) <0.0001)   goto L185;
              TO=T;
            }

          goto fin;


L185: if (T<13.){ if (2.*EE*(1.-2.*EPS*T-V*(3.-12.*T*T+4.*pow(T,4))*exp(-T*T)) >0.0) goto L205; }
          else      { if (2.*EE*(1.-2.*EPS*T) >0.0) goto L205; }


L190:

          for(I=1; I<=20; I++) {
              TO=double(I)/10.;
              GL=EE*(1.-EPS*TO-V*(3.-2.*TO*TO)*exp(-TO*TO));
              if (GL>=0.0) goto L175;
            }

L205:
          SHLL=EE*pow(T,2)-FF*pow(T,3)+SSHELL*(1.-2.*pow(T,2))*exp(-pow(T,2));
        }
      //-----------------

      YMASS=SMASS+SHLL*(1.-_NOSHL);

      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
    }  else {
      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

      if(A== 1.&&Z== 1.)YMASS=7.29;
      else if(A== 1.&&Z== 0.)YMASS=8.07;
      else if(A== 2.&&Z== 1.)YMASS=13.14;
      else if(A== 3.&&Z== 1.)YMASS=14.95;
      else if(A== 3.&&Z== 2.)YMASS=14.93;
      else if(A== 4.&&Z== 2.)YMASS=2.43;
      else if(A== 5.&&Z== 2.)YMASS=11.39;
      else if(A== 5.&&Z== 3.)YMASS=11.68;
      else if(A== 6.&&Z== 3.)YMASS=14.09;
      else if(A== 7.&&Z== 3.)YMASS=14.91;
      else if(A== 8.&&Z== 3.)YMASS=20.95;
      else if(A== 9.&&Z== 3.)YMASS=24.95;
      else if(A== 6.&&Z== 4.)YMASS=18.37;
      else if(A== 7.&&Z== 4.)YMASS=15.77;
      else if(A== 8.&&Z== 4.)YMASS=4.94;
      else if(A== 9.&&Z== 4.)YMASS=11.35;
      else if(A==10.&&Z== 4.)YMASS=12.61;
      else if(A==11.&&Z== 4.)YMASS=20.17;
      else if(A== 7.&&Z== 5.)YMASS=27.94;
      else if(A== 8.&&Z== 5.)YMASS=22.92;
      else if(A== 9.&&Z== 5.)YMASS=12.42;
      else if(A==10.&&Z== 5.)YMASS=12.05;
      else if(A==11.&&Z== 5.)YMASS=8.67;
      else if(A==12.&&Z== 5.)YMASS=13.37;
      else if(A==13.&&Z== 5.)YMASS=16.56;
      else if(A== 9.&&Z== 6.)YMASS=28.91;
      else if(A==10.&&Z== 6.)YMASS=15.70;
      else if(A==11.&&Z== 6.)YMASS=10.65;
      else if(A==12.&&Z== 6.)YMASS=0.;
      else if(A==13.&&Z== 6.)YMASS=3.13;
      else if(A==14.&&Z== 6.)YMASS=3.02;
      else if(A==15.&&Z== 6.)YMASS=9.87;
      else if(A==16.&&Z== 6.)YMASS=13.69;
      else if(A==12.&&Z== 7.)YMASS=17.34;
      else if(A==13.&&Z== 7.)YMASS=5.34;
      else if(A==14.&&Z== 7.)YMASS=2.86;
      else if(A==15.&&Z== 7.)YMASS=0.10;
      else if(A==16.&&Z== 7.)YMASS=5.68;
      else if(A==17.&&Z== 7.)YMASS=7.87;
      else if(A==18.&&Z== 7.)YMASS=13.27;
      else if(A==13.&&Z== 8.)YMASS=23.10;
      else if(A==14.&&Z== 8.)YMASS=8.01;
      else if(A==15.&&Z== 8.)YMASS=2.85;
      else if(A==16.&&Z== 8.)YMASS=-4.74;
      else if(A==17.&&Z== 8.)YMASS=-.81;
      else if(A==18.&&Z== 8.)YMASS=-.78;
      else if(A==20.&&Z== 8.)YMASS=3.80;
      else if(A==21.&&Z== 8.)YMASS=8.12;
      else if(A==15.&&Z== 9.)YMASS=17.66;
      else if(A==16.&&Z== 9.)YMASS=10.69;
      else if(A==17.&&Z== 9.)YMASS=1.95;
      else if(A==18.&&Z== 9.)YMASS=0.87;
      else if(A==19.&&Z== 9.)YMASS=-3.43;
      else if(A==20.&&Z== 9.)YMASS=-0.017;
      else if(A==21.&&Z== 9.)YMASS=-0.047;
      else if(A==22.&&Z== 9.)YMASS=2.83;
      else if(A==23.&&Z== 9.)YMASS=3.35;
      else if(A==20.&&Z==10.)YMASS=-7.04;
      else if(A==21.&&Z==10.)YMASS=-5.73;
      else if(A==22.&&Z==10.)YMASS=-8.03;
      else if(A==23.&&Z==10.)YMASS=-5.15;
      else if(A==24.&&Z==10.)YMASS=-5.95;
      else if(A==25.&&Z==10.)YMASS=-2.15;
      else if(A==26.&&Z==10.)YMASS=-0.019;
      else if(A==20.&&Z==11.)YMASS=6.84;
      else if(A==21.&&Z==11.)YMASS=-2.18;
      else if(A==22.&&Z==11.)YMASS=-5.18;
      else if(A==23.&&Z==11.)YMASS=-9.53;
      else if(A==24.&&Z==11.)YMASS=-8.42;
      else if(A==25.&&Z==11.)YMASS=-9.36;
      else if(A==26.&&Z==11.)YMASS=-6.89;
      else if(A==27.&&Z==11.)YMASS=-5.63;
      else if(A==22.&&Z==12.)YMASS=-0.39;
      else if(A==23.&&Z==12.)YMASS=-5.47;
      else if(A==24.&&Z==12.)YMASS=-13.93;
      else if(A==25.&&Z==12.)YMASS=-13.19;
      else if(A==26.&&Z==12.)YMASS=-16.21;
      else if(A==27.&&Z==12.)YMASS=-14.58;
      else if(A==28.&&Z==12.)YMASS=-15.02;
      else if(A==29.&&Z==12.)YMASS=-10.75;
      else if(A==24.&&Z==13.)YMASS=-0.052;
      else if(A==25.&&Z==13.)YMASS=-8.91;
      else if(A==26.&&Z==13.)YMASS=-12.21;
      else if(A==27.&&Z==13.)YMASS=-17.2;
      else if(A==28.&&Z==13.)YMASS=-16.85;
      else if(A==29.&&Z==13.)YMASS=-18.21;
      else if(A==30.&&Z==13.)YMASS=-15.89;
      else if(A==31.&&Z==13.)YMASS=-15.10;
      else if(A==26.&&Z==14.)YMASS=-7.14;
      else if(A==27.&&Z==14.)YMASS=-12.38;
      else if(A==28.&&Z==14.)YMASS=-21.49;
      else if(A==29.&&Z==14.)YMASS=-21.89;
      else if(A==30.&&Z==14.)YMASS=-24.44;
      else if(A==31.&&Z==14.)YMASS=-22.95;
      else if(A==32.&&Z==14.)YMASS=-24.09;
      else if(A==33.&&Z==14.)YMASS=-20.57;
      else if(A==28.&&Z==15.)YMASS=-7.16;
      else if(A==29.&&Z==15.)YMASS=-16.95;
      else if(A==30.&&Z==15.)YMASS=-20.20;
      else if(A==31.&&Z==15.)YMASS=-24.44;
      else if(A==32.&&Z==15.)YMASS=-24.30;
      else if(A==33.&&Z==15.)YMASS=-24.94;
      else if(A==34.&&Z==15.)YMASS=-24.55;
      else if(A==35.&&Z==15.)YMASS=-24.94;
      else if(A==30.&&Z==16.)YMASS=-14.06;
      else if(A==31.&&Z==16.)YMASS=-19.04;
      else if(A==32.&&Z==16.)YMASS=-26.01;
      else if(A==33.&&Z==16.)YMASS=-26.58;
      else if(A==34.&&Z==16.)YMASS=-29.93;
      else if(A==35.&&Z==16.)YMASS=-28.85;
      else if(A==36.&&Z==16.)YMASS=-30.66;
      else if(A==37.&&Z==16.)YMASS=-26.91;
      else if(A==38.&&Z==16.)YMASS=-26.86;
      else if(A==32.&&Z==17.)YMASS=-13.33;
      else if(A==33.&&Z==17.)YMASS=-21.00;
      else if(A==34.&&Z==17.)YMASS=-24.44;
      else if(A==35.&&Z==17.)YMASS=-29.01;
      else if(A==36.&&Z==17.)YMASS=-29.52;
      else if(A==37.&&Z==17.)YMASS=-31.77;
      else if(A==38.&&Z==17.)YMASS=-29.80;
      else if(A==39.&&Z==17.)YMASS=-29.80;

      if(A==40.&&Z==17.)YMASS=-27.54;
      else if(A==34.&&Z==18.)YMASS=-18.38;
      else if(A==35.&&Z==18.)YMASS=-23.05;
      else if(A==36.&&Z==18.)YMASS=-30.23;
      else if(A==37.&&Z==18.)YMASS=-30.95;
      else if(A==38.&&Z==18.)YMASS=-34.72;
      else if(A==39.&&Z==18.)YMASS=-33.24;
      else if(A==40.&&Z==18.)YMASS=-35.04;
      else if(A==41.&&Z==18.)YMASS=-33.07;
      else if(A==42.&&Z==18.)YMASS=-34.42;
      else if(A==43.&&Z==18.)YMASS=-31.98;
      else if(A==44.&&Z==18.)YMASS=-32.27;
      else if(A==36.&&Z==19.)YMASS=-17.43;
      else if(A==37.&&Z==19.)YMASS=-24.80;
      else if(A==38.&&Z==19.)YMASS=-28.80;
      else if(A==39.&&Z==19.)YMASS=-33.80;
      else if(A==40.&&Z==19.)YMASS=-33.53;
      else if(A==41.&&Z==19.)YMASS=-35.55;
      else if(A==42.&&Z==19.)YMASS=-35.02;
      else if(A==43.&&Z==19.)YMASS=-36.59;
      else if(A==44.&&Z==19.)YMASS=-35.81;
      else if(A==45.&&Z==19.)YMASS=-36.61;
      else if(A==46.&&Z==19.)YMASS=-35.42;
      else if(A==47.&&Z==19.)YMASS=-35.69;
      else if(A==48.&&Z==19.)YMASS=-32.22;
      else if(A==39.&&Z==20.)YMASS=-27.28;
      else if(A==40.&&Z==20.)YMASS=-34.85;
      else if(A==41.&&Z==20.)YMASS=-35.14;
      else if(A==42.&&Z==20.)YMASS=-38.54;
      else if(A==43.&&Z==20.)YMASS=-38.40;
      else if(A==44.&&Z==20.)YMASS=-41.46;
      else if(A==45.&&Z==20.)YMASS=-40.81;
      else if(A==46.&&Z==20.)YMASS=-43.14;
      else if(A==47.&&Z==20.)YMASS=-42.34;
      else if(A==48.&&Z==20.)YMASS=-44.22;
      else if(A==49.&&Z==20.)YMASS=-41.28;
      else if(A==50.&&Z==20.)YMASS=-39.57;
    }

fin:
  //if(YMASS==999) fprintf(f09,"\n A=%.0f Z=%.0f YMASS=999 !",A,Z);

  return YMASS;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

