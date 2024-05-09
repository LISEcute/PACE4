#include "ParticleStates.h"
//#include <time.h>
#include "fisrot.h"
#include <iostream>
#include <QDebug>
#include <QApplication>
#include <QProgressDialog>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tfusion.h"
#include "ftype.h"
#include <QTextStream>
#include <QDate>
#include <QTime>


#include "pace.h"
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))

#define HIGH_WORD(k) ((short int)(k >> 16))
#define LOW_WORD(k)  ((short int)(k & 0xFFFF))

///   it is necessary to clear each time  CHECK before run PACE3 !!!!!!
///-------------------- static to clear
int _cl_K,_cl_IPRNT,_cl_KK,_cl_IZPREV,_cl_INPREV;
//---------------------


// OLEG
//---only for READ
double _QCNr, _EXPSIGr, _BARFACr, _AGRAZr;
double  _FGE1r,_FGM1r,_FGE2r,_FGM2r;
double _LowLimit, _HighLimit;
int ParticleFlag;
//---NEW one from 27-01-03

int _BarrierMode=0; //  = 0 classic
//   = 1 quantum-mechanic
double _h_omega = 3;  // 3 MeV - default value

int _BatchMode = 0; // 0-No, 1-Yes
int _CurrentBatchStep=0;
int _BatchNsteps=10;

//---also for WORK

int _INPUT, _IDIST;
int _IZP,_IAP,_IZT,_IAT, _IZC, _IAC;
double _SPr, _ST,  _QCN;
double _ELAB,_ExpSig,_AGRAZ,_ELOSS, _ELAB_MAX ;

int _LMINN, _JCMAX,_JCMAXr;
int   _NPOT;   // number of potentials
//int _ITAKE[5];
double _EEXCN,_EEXCN_MAX, _EREC, _AJNUC;
//double _FISSBR;
//=========      common/MASS/
double _FLMASS,_BARFAC,_ARATIO,_BARMAX,_FBARR[111];

int _IZR[5],_INR[5], _LFISS;

int _IZPART[5]= {0,0,1,2,0};
int _INPART[5]=  {0,1,0,2,0};

//=========      common/QEJS/
int _MEBIN[5];

//=========      common/GAM/
double _FGE1,_FGM1,_FGE2,_FGM2,_DISCPR,_ERGMIN ,_EDL[37][7];
int  _IADEF,_NDISC,_NDL[7],_IZDL[7],_INDL[7],_JDL[37][7];

//=========      common/CP/
double _ALEV[Max_MOM+1];
extern double _ALEVB[Max_MOM+1];

//=========      common/SRCH/
double _PROB[5556],_EMAX[5],_EMIN[5],_BE[5],_GE1,_GE2,_GM1,_GM2;
int    _IPROB,_MJ[5556],_IMODL;

//=========      common/POT/
double _V0D[5],_R0RD[5],_ARD[5],_R0CD[5],_RCCC[4],_VQS[4],_VQ0[4],_W0D[5],_R0ID[5],_AID[5],
_RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
int _NPD[5],_IMAG[5],_IRAD[5];

//=========      common/SCH/
double _EBIN[5][DIM_EBIN];//,_RMASS[5];  //MPK    // OLEG
double _RMASS[5] = {0};
int    _MAXJ[5][DIM_EBIN],_MAXJS[5][DIM_EBIN]; // OLEG

//=========      common/DAT/
double _SPART[5]={0,0.5,.5,0.,0.};
double _PMASS[5]={0,8.07169,7.28922,2.42494,0.};

//=========      common/OUT/


double _RLEV[4][DIM_RLEV];        //   base 0!

//=========      common /XQSIG/
double _SIGMA;
int _NCASC;
double _DSIG[19];
//=========      common/QMMJ/
double _FYRST,_FACLA,_EROT[Max_MOM+1];
int    _MAXC=0;

//=========      common/XQRES/
double _ENERGY,_VZC;
int    _IZ,_IN,_ITRAC;

//=========      common/XQANG/
double _SUM[11][301],_CT;     // SUM - transposed
int    _MDIR;

//=========      common/XQSHL/
int _NOSHL=0;

//=========      common/QSIER/
double _BARLD[81];
int  _IS;

//=========      common/FAC/
double FA[204];
//==================     common/EM/
int eNT[5],eN[4][33][19],eNG[33],eNW[4][5][19],eNGW[5],eNGB[33],eNGBW[5],eNGT[33],eNGTW[33];
double eEW[4],eEW1[5],eEW2[4];

//     ***********
//      MONITORING OF OUTGOING PARTICLES.
//      common/MOM/

//==================     common /XQTT/
double _UX,_EX,_TT;


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw


//  all arrays for base 1 !!!

FILE *f10, *f09, *f02, *f05, *f_particles;

QTextStream *htmlStream;
QFile *htmlFile;

//extern char *s35x;
const char s35x[10] = "         ";
extern char *s100m;
//char *s100m = "      ";
const char *s40star="************************************";
const char *s10star="***********";



extern double RANF();
extern double MASSES(int IZ, int IN, int &opt);
extern void PLM();
extern void TLLL(double ER, int NP,int IZ,int IN, double TL[25], int &LLMAX, QTextStream& s);
extern void FIND(double E,int &IXRSM, double *EBIN, int MEBIN, int MODE, QTextStream& s);
extern void SEARCH(fusion_event *a, fusion_event &c, int &MODE, QTextStream& m);
extern void MJRAN(fusion_event &csi, fusion_event &csf, int MODE, double EP, QTextStream& s);
extern void PRODCT(int II,  int IZ, int IN, int NCASC, fusion_event **a, const char *filename_cs, QTextStream& s);
extern void STATIS(fusion_event **a, const char *filename_evt, ClassParticleFlags &ParticleFlag, QTextStream& s);


/**** MICHELLE ADD FROM tfusion2a2.cpp *****************/
double _ALEVB[Max_MOM+1] = {0};


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//int tfusion::PACE3(const char *filename_rtf, const char *filename_evt, const char *filename_cs, const char *filename_particles)//, QString html_results)
int PACE::PACE3(const char *filename_rtf, const char *filename_evt,
                const char *filename_cs, const char *filename_particles,
                const char *filename_html)
{

  short int IBUF2[10];
  //MK    QVector<int> buf2List(10);
  //MK    QVector<QVector<int> > buf2Link;

  float BUF3[5];
  //MK    QVector<float> buf3List(10);
  //MK    QVector<QVector<float> > buf3Link;

  double  PROBM[5], EXRSM[5], TL[25];
  int     IXRSM[5], IXPR[5], IXMIN[5], NSPC[5][52] ;     //NOT transported
  int I, L, K, J, IST, II;
  double SLEV, Fprob=0;   // FPROB - Fission probability

  double  ERGNOT=0, EP=0, AC, AP, FISSBR,EU,SC;
  int MODE=2;
  int IZMAX,INMAX,IA,LA,MA,MERG,LERG,IDUM,IKQ,JKQ,IAC;
  double   A3,RD,VCOULP,VCOULA,DISC,t;
  double A, Z, AL, AAA,AJN, ETRAN=0,EDUM;
  // extern bool INTRAC;
  char  NTITL_B[2][5]={"VOL ","SURF"};
  double SPROB=0,TPROB=0,ECLOSS=0, EXMIN=0;
  //fusion_event  csi, c, csf, **a;
  fusion_event  csi, csf, **a;
  fusion_event c;

  //     time_t ti = time(0);      v.4.34.2
  //    struct tm *mytime = localtime(&ti);
  QDate mydate = QDate::currentDate();
  QTime mytime = QTime::currentTime();


  double Emin_wait;
  double CoefBarfac;
  double SpinBarFac;


  if(_NCASC>Max_NCASC) _NCASC=Max_NCASC;

  a = new fusion_event*[_NCASC+1];
  for(I=0; I<=_NCASC; I++) a[I]=new fusion_event;

  _QCN=_QCNr;
  _FGE1 = _FGE1r;
  _FGM1 = _FGM1r;
  _FGE2 = _FGE2r;
  _FGM2 = _FGM2r;
  _BARFAC = _BARFACr;
  _ExpSig = _EXPSIGr;
  _AGRAZ  = _AGRAZr;
  _JCMAX  = _JCMAXr;
  for(I=0; I<=Max_MOM; I++) _ALEV[I]=_ALEVB[I];

  _cl_K=_cl_IPRNT=_cl_KK =_cl_IZPREV=_cl_INPREV=0;

  f09=fopen (filename_rtf,"wt");

  f02=f_particles=nullptr;

  QString html_filename(filename_html);
  QFile *htmlFile = new QFile(html_filename);
  if (!htmlFile->open(QIODevice::WriteOnly)){
      qDebug() << "COULD NOT OPEN HTML FILE!"  << html_filename;
      return -1;
    }

  htmlStream = new QTextStream(htmlFile);//,QIODevice::ReadWrite);

  if(_IDIST!=0)
    {
      f02=fopen(filename_evt,"wb");
    }

  if(ParticleFlags.CanBeOpenedParticleFile())
    {
      f_particles = fopen (filename_particles,"wt");
      ParticleFlags.SetParticleState(ep_opened);
    }

  _BARMAX=30.;

  //==========================================================================
  //     NO FISSION CALCULATION FOR B.F. ABOVE BARMAX

  PLM();
  int MODES=4;
  _NPOT=(_INPUT==1 ? 4 : 3);   // number of potentials
  _IMODL=4*DIM_RLEV;
  _BE[4]=0.;
  ECLOSS=0.;
  _EMIN[1]=.01;
  //     DISCPR IS DISCRIMINATOR ON PROBABILITY VECTOR (PROB) TO PREVENT
  //            SEARCH FOR HIGHLY UNLIKELY EVENTS. CHANGE TO 0.001 IF YOU
  //            WANT RARE EVENTS AND ARE USING GOOD STATISTICS.
  _EMIN[4]=0.;

  fprintf(f09,
          "{\\rtf1\\ansi\\deff0\\deftab708{\\fonttbl{\\f0\\fmodern\\fprq1 Courier New;}{\\f1\\froman Symbol;}{\\f2\\froman Times New Roman;}}"
          "{\\colortbl ;\\red255\\green0\\blue0;\\red51\\green51\\blue153;\\red0\\green0\\blue255;\\red60\\green170\\blue140;\\red128\\green128\\blue127;}"
          "\\viewkind4\\uc1\\pard\\tx450\\tx1800");


  fprintf(f09,"\\f0\\fs10\\cf0\\par%s \\par ",filename_rtf);

  QString qbuf = QString::asprintf("%s ** %02d:%02d ** %02d-%02d-%04d",
                                   PACE_version, mytime.hour(), mytime.minute(),
                                   mydate.day(), mydate.month(), mydate.year());


  *htmlStream << "<!DOCTYPE html><html><head><title>PACE4 Results</title><style>td{padding:2px;} h3#mode{color:blue; margin-left:650px; }</style></head><body>";
  *htmlStream << "<p>" + html_filename + "<BR>" + QString(filename_rtf) + "</p>";
  *htmlStream << "<p>" + qbuf + "</p>";

  fprintf(f09,"%s \\par\\cf0\\fs16 ",qbuf.toStdString().c_str());

  fprintf(f09,"\\f0\\par"
              "\\fs22%s%s%s\\fs40\\b\\cf1\\f2  P A C E   4\\par"
              "\\f0\\par"
              "\\fs22%s%s%s\\cf4 modified JULIAN \\fs16\\par"
              "\\par\\f0\\cf0\\b0\\fs16"
              "\\cf3%s projection angular-momentum coupled evaporation Monte Carlo code  %s\\par"
              "%s angular distributions obtained using M-states of angular momentum %s\\par"
              "\\cf0\\par"
          ,
          s35x,s35x,s35x,
          s35x,
          s35x,s35x,
          s10star,s10star,s10star,s10star
          );
  *htmlStream << "<div align=center><h1 style=\"color:red\">P A C E 4</h1><h3 style=\"color:green\">modified JULIAN</h3><h3 style=\"color:blue\"> "
                 "projection angular-momentum coupled evaporation Monte Carlo code  </h3><h3 style=\"color:blue\">"
                 " angular distributions obtained using M-states of angular momentum</h3></div><h3 id=\"mode\">";
  fprintf(f09,"\\b\\cf3\\fs20 %s%s%s%s%s%s%s\\f2 MODE=%d\\f0\\cf0\\b0\\fs16",s35x,s35x,s35x,s35x,s35x,s35x,s35x,_INPUT);
  *htmlStream  << "MODE = " << QString::number(_INPUT) << "</h3>";

  int     IBOUND=0;
  int     INOTGS=0;
  int     INOTG1=0;
  double  SNOTGS=0.;
  double  ENOTGS=0.;
  int     NODES=MODES-1;

  //      ARATIO IS (LITTLE-A-SADDLE)/(LITTLE-A-G.S.)
  //      DEFAULT FOR NO INPUT IS 1.000

  if(_ARATIO==0) _ARATIO=1;
  double ACASC=_NCASC;
  ACASC=sqrt(ACASC)*3.0+10.;

  if(ACASC<500.)ACASC=500.;
  _DISCPR=1./ACASC;

  int countMax = 25; //initialize now, set later
  int progress_inc = 0;
  int val=0;

  QProgressDialog progress("Calculating..","Cancel",0,countMax);
  progress.setWindowFlags(Qt::CustomizeWindowHint |
                          Qt::WindowCloseButtonHint);
  progress.setWindowIcon(QIcon(":/pace4.png"));
  progress.setWindowModality(Qt::WindowModal);
  progress.show();

  if( COMPOS(_IZ,_IN,_ENERGY,_VZC,_MAXC,_SIGMA,ECLOSS,_INPUT,*htmlStream)<0)
    {
      if(_IDIST!=0 && f02)  {fclose(f02); f02=nullptr;}
      if(f_particles)
        {
          fclose(f_particles);
          f_particles=nullptr;
          //ParticleFlags.ClearParticleState(ep_opened);
        }
      goto BadEnd;
    }

  //      IADEF=2 IS ASSUMED SPHERICAL NUCLEI FOR LEVEL DENSITIES.
  //      IADEF=1 IS CONSIDERED A DEFORMED NUCLEUS.
  //      THE DECISION IS MADE ONLY FOR THE Z,N OF THE INITIAL NUCLEUS.
  //     COMPOS FORMS INITIAL SPIN DISTRIBUTION DEPENDING ON INPUT

  fprintf(f09,"\\par\\cf2\\f2\\fs20  *** ");   *htmlStream << "<p> ***";
  if( (_IN>=77 && _IN<=110) || _IN>138)
    {
      _IADEF=1;
      fprintf(f09,"Deformed nucleus. Backbending at J=20 hbar");
      *htmlStream << "Deformed nucleus. Backbending at J=20 hbar";
    }
  else {
      _IADEF=2;
      fprintf(f09,"Spherical nucleus level density");
      *htmlStream << "Spherical nucleus level density";
    }

  fprintf(f09,"\\cf0\\f0\\fs16");   *htmlStream << "</p>";

  IZMAX=_IZ;
  INMAX=_IN;
  int opt1;
  _FLMASS = MASSES(_IZ,_IN, opt1);
  Z=IZMAX;


            if(Z==68)_IS=1;
  else      if(Z==76)_IS=2;
  else      if(Z==84)_IS=3;
  else               _IS=0;

  AAA=IZMAX+INMAX;
  fprintf(f09,"\\cf2\\f2\\fs20 ");

       if(_IS==1)AAA=158.;
  else if(_IS==2)AAA=186.;
  else if(_IS==3)AAA=210.;
//  else if(_IS==0){
//OT 07/26/23      fprintf(f09,"\\par  *** Sierk barrier not found for A =%5.0f  Z=%4.0f",AAA,Z);
//OT 07/26/23      *htmlStream << "<p>  *** Sierk barrier not found for A =" << QString::number(AAA) << " Z= " << QString::number(Z) << "</p>";
//    }

  //MM   IF(BARFAC.LT.(-0.01))IS=6

  if(_BARFAC==0) _BARFAC=-1;

  CoefBarfac = (_BARFAC < 0 ? fabs(_BARFAC) : 1.);
  SpinBarFac = (_BARFAC > 0 ? _BARFAC : 0);

  _BARFAC=fabs(_BARFAC);

  //     BARFAC IS FACTOR ONLY
  //      if(_IS>0)    OT
  if(CoefBarfac!=1)
    {
      fprintf(f09,"\\par  *** Sierk barriers X factor of %-9.2f",CoefBarfac);
      *htmlStream << "<p>  *** Sierk barriers X factor of " << QString::number(CoefBarfac) << "</p>";
    }

  //     BARFAC MULTIPLIES THE C.P.S. BARRIER TO GET YOUR OWN.
  for(II=1; II<=80; II++)
    {
      AL=II-1;
      _BARLD[II]= FISROT(AAA,Z,AL,1.);
    }

  A=IZMAX+INMAX;

  FISSBR =  SIERK(A,Z,SpinBarFac,CoefBarfac);

  fprintf(f09,"\\par\\par  *** Input fission barrier = %.2f MeV at L=0  taken from Sierk"
              "\\par  *** G.S. little A multiplied by factor %.3f obtain saddle level density"
              "\\par  *** No fission calculation for barrier above %.1f MeV",
          FISSBR,_ARATIO,_BARMAX);

  *htmlStream << "<p>&nbsp;</p><p> *** Input fission barrier = " << QString::number(FISSBR,'g',4) << " MeV at L=0 taken from Sierk</p>" <<
                 "<p> *** G.S. little A multiplied by factor " << QString::number(_ARATIO,'g',4) << " to obtain saddle level density</p>" <<
                 "<p> *** No fission calculation for barrier above " << QString::number(_BARMAX,'g',3) << " MeV</p>";
  if(FISSBR==0)
    {
      fprintf(f09,"\\par  *** G.S. level density little a multiplied by %.3f to "
                  "obtain saddle point level density for fission calculation"
                  "\\par  *** Zero spin fission barrier is %6.2f MeV",_ARATIO,FISSBR);
      *htmlStream << "<p> *** G.S. level density little a multiplied by " << QString::number(_ARATIO,'g',4) << " to "
                                                                                                               "obtain saddle point level density for fission calculation</p>"
                                                                                                               "<p> *** Zero spin fission barrier is " << QString::number(FISSBR,'g',3) << "</p>";
    }

  if(_FACLA==0.)
    {
      fprintf(f09,"\\par  *** FACLA=0 ----  Gilbert-Cameron little-A selected\\par");
      *htmlStream << "<p> *** FACLA=0 ----  Gilbert-Cameron little-A selected </p>";
    }
  else if(_FACLA >0.)
    {
      fprintf(f09,"\\par  *** Little-A = MASS / %4.1f\\par",_FACLA);
      *htmlStream << "<p>   *** Little-A = MASS / " << QString::number(_FACLA,'g',3) << "</p>";
    }


  if(_NOSHL==1)
    {
      fprintf(f09,"\\par  *** Lysekil masses with no shell correction used \\par%s\\par", s40star);
      *htmlStream << "<p> *** Lysekil masses with no shell correction used</p>";
    }
  fprintf(f09,"\\cf0\\f0\\fs16");

  A3=pow(A,.3333);

  RD=A3+1.;
  VCOULP=double(IZMAX)/RD;

  RD=A3+1.587;
  VCOULA=2.*double(IZMAX)/RD;
  DISC  =1.5+6.5*exp(-0.01*A);
  //     ***************************
  _EMAX[1]=40.;
  _EMAX[2]=VCOULP*DISC;
  _EMAX[3]=VCOULA*DISC;
  _EMAX[4]=20.;

  _EMIN[2]=0.015*Z*(1.+0.01*Z);
  _EMIN[3]=2.*_EMIN[2];
  //      EMIN ARE THE SMALLEST ENERGIES FOR WHICH RELIABLE CROSS SECTIONS
  //          ARE CALCULABLE.  THE FORMULAE WERE OBTAINED FROM EXPERIENCE.
  //     EMAX ARE MAXIMUM EMISSION ENERGIES SO AS NOT TO WASTE TIME
  //     COMPUTING UNLIKELY EVENTS
  for(IKQ=1; IKQ<=4; IKQ++)
    for(JKQ=1; JKQ<=51; JKQ++)
      NSPC[IKQ][JKQ]=0;

  fprintf(f09,"\\par\\par"
              "\\par\\b      Energy range for      neutron    proton      alpha      gamma\\b0");
  *htmlStream << "<p>&nbsp;</p><table cellpadding='5'><tr><th width=15% align=left>Energy range for</th><th>neutron</th><th>proton</th><th>alpha</th><th>gamma</th></tr>";
  fprintf(f09,"\\par      minimal               %6.2f     %6.2f     %6.2f     %6.2f",
          _EMIN[1],_EMIN[2],_EMIN[3],_EMIN[4]);
  *htmlStream << "<tr><td>minimal</td><td>" << QString::number(_EMIN[1],'f',2) << "</td><td>" << QString::number(_EMIN[2],'f',2) << "</td><td>"
                                                                                                                                 <<  QString::number(_EMIN[3],'f',2) << "</td><td>" << QString::number(_EMIN[4],'f',2)
      << "</td></tr>";
  fprintf(f09,"\\par      minimal               %6.2f     %6.2f     %6.2f     %6.2f",
          _EMAX[1],_EMAX[2],_EMAX[3],_EMAX[4]);
  *htmlStream << "<tr><td>minimal</td><td>" << QString::number(_EMAX[1],'f',2) << "</td><td>" <<  QString::number(_EMAX[2],'f',2) << "</td><td>"
                                                                                                                                  << QString::number(_EMAX[3],'f',2) << "</td><td>" << QString::number(_EMAX[4],'f',2) << "</td></tr></table>";
  fprintf(f09,"\\par\\par\\cf5\\f2\\fs20  *** Internal probability discriminator of program set to %.4f \\cf1    Number of cascades is %d\\cf0\\f0\\fs16",
          _DISCPR,_NCASC);
  *htmlStream << "<p>&nbsp;</p><p style=\"color:#555\">*** Internal probability discriminator of program set to " << QString::number(_DISCPR)
              << "</p><h4 style=\"color:red\">Number of cascades is " << QString::number(_NCASC) << "</h4>";
  //     OUTEM IS PARTICLE EMISSION OUTPUT SUBROUTINE
  htmlStream->flush();
  if(_IDIST!=0)  OUTEM(1,*htmlStream);      // initialisation

  fprintf(f09,"\\par\\par\\b       Optical model parameters for light emitted particles\\b0"
              "\\par\\b\\fs14    V       *E     *E**2    R0R    ARD    R0C      W0     *E  "
              "*E**2    R01    AID  RMCHD   NPD  IMAG IRAD\\b0");
  *htmlStream << "<p>&nbsp;</p><table cellspacing=3><tr><th colspan=15>Optical model parameters for light emitted particles</th></tr>"
              << "<tr><th>V</th><th>*E</th><th>*E**2</th><th>R0R</th><th>ARD</th><th>R0C</th><th>W0</th><th>*E</th><th>*E**2</th><th>R01</th>"
                 "<th>AID</th><th>RMCHD</th><th>NPD</th><th>IMAG</th><th>IRAD</th></tr>";
  //     ***************************
  for(I=1; I<=3; I++) {
      if(!_ITAKE[I] || _V0D[I]==0) OPTPOT(I,IZMAX,INMAX);
      fprintf(f09,"\\par %7.3f %7.3f %7.3f %6.3f %6.3f %6.3f %7.3f %6.3f %6.3f %6.3f %6.3f %6.3f %5d  %s %3d",
              _V0D[I],_V01D[I],_V02D[I],_R0RD[I],_ARD[I],_R0CD[I],
              _W0D[I],_W01D[I],_W02D[I],_R0ID[I],_AID[I],_RMCHD[I],
              _NPD[I],NTITL_B[_IMAG[I]],_IRAD[I]);
      *htmlStream << "<tr><td>" << QString::number(_V0D[I],'f',3) << "</td><td>" << QString::number(_V01D[I],'f',3) <<"</td><td>" << QString::number(_V02D[I],'f',3) << "</td><td>" <<QString::number(_R0RD[I],'f',3)
                  <<"</td><td>" << QString::number(_ARD[I],'f',3) <<"</td><td>" << QString::number(_R0CD[I],'f',3) <<"</td><td>" << QString::number(_W0D[I],'f',3) <<"</td><td>" << QString::number(_W01D[I],'f',3)
                 << "</td><td>" <<QString::number(_W02D[I],'f',3)<<"</td><td>" << QString::number(_R0ID[I]) << "</td><td>" << QString::number(_AID[I]) <<"</td><td>" << QString::number(_RMCHD[I],'f',3)
                 <<"</td><td>" << QString::number(_NPD[I],'f',3) << "</td><td>" <<QString::fromUtf8(NTITL_B[_IMAG[I]]) << "</td><td>" <<QString::number(_IRAD[I],'f',3) << "</td></tr>";
    }

  *htmlStream << "</table>";
  fprintf(f09,"\\fs16" );
  EDUM=0;     IDUM=0;
  TLLL(EDUM,IDUM,IZMAX,INMAX,TL,IDUM, *htmlStream);

  //     TLLL PREPARES ALL N, P, AND ALPHA TRANSMISSION COEFFICIENTS IN ADVANCE.
  //     THE CHANGE IN TRANSMISSION COEFFICIENTS FOR SECOND AND SUBSEQUENT DECAYS
  //     IS SIMULATED BY SHIFTING THE ENERGY AT WHICH
  //     THE TL-S ARE CALCULATED ACCORDING TO THE SHIFT OF THE COULOMB
  //     BARRIER OF THE SECONDARY RESIDUAL NUCLEI.
  //     TL-S AT THIS ENERGY ARE CALCULATED BY INTERPOLATION
  //
  //     **************************

  //L46:  fscanf (10,41)FGE1,FGM1,FGE2,FGM2;

  if(_FGE1==0 ){   // OLEG  ADD  ALSO CHECK CARD
      if(A<=100) { _FGE1=.00008;      _FGE2=4.8;      _FGM1=.025;    _FGM2=.0195;}
      else if(A<=126) { _FGE1=.000014;     _FGE2=5.9;      _FGM1=.01;     _FGM2=.00088;}
      else if(A<=150) { _FGE1=.000046;     _FGE2=7.7;      _FGM1=.007;    _FGM2=.058;}
      else            { _FGE1=.000011;     _FGE2=9.0;      _FGM1=.010;    _FGM2=1.2;}
    }
  fprintf(f09,"\\par\\par\\b         E.M.Transition strengths in Weisskopf units\\b0\\fs14"
              "\\par        E1 = %-10f   M1 = %-10f   E2 = %-10f    M2 =  %-10f\\fs16",
          _FGE1, _FGM1, _FGE2, _FGM2);
  *htmlStream << "<p>&nbsp;</p><table cellspacing=10><tr><th colspan=4>E.M.Transition strengths in Weisskopf units</th></tr>"
                 "<tr><td> E1 = " << QString::number(_FGE1,'f',6) << "</td><td>   M1 = " << QString::number(_FGM1,'f',6)
              << "</td><td>   E2 = " << QString::number(_FGE2,'f',6) << "</td><td>   M2 = " << QString::number(_FGM2,'f',6) << "</td></tr></table>";
  _NDISC=1;


  //      DISCRETE LEVELS (IF DESIRED) FOR VARIOUS LEVEL DENSITIES
  //      ARE NOW READ IN. THE NUMBER OF FINAL NUCLEI
  //      THAT HAVE DISCRETE LEVELS IS NDISC.
  //      NOTE *****  NDISC SHOULD NOT BE MORE THAN 6.
  //           *****  THE MAXIMUM NUMBER OF LEVELS FOR EACH NUCLEUS IS 36
  /*
do {
//     fscanf (10,501,IOSTAT=IEOF)NDL[_NDISC],IZDL[_NDISC],INDL[_NDISC];
      if(IEOF!=0)           { NDISC--;  break;}
      NLL=NDL[_NDISC];
      if( NLL<=0 || NLL>36) { NDISC--;  break;}

//      fscanf (10,504)(EDL[I][NDISC],JDL[I][NDISC],I=1,NLL);
      fprintf(f09,"\\par  Discrete levets.  Z N =%d %d  develop print LATER", IZDL[_NDISC],_INDL[NDISC]);
      //',2I4/(8(F6.2,I4,4X)/));)IZDL[NDISC],INDL[NDISC],(EDL[I][NDISC],JDL[I][NDISC],I=1,N LL);
      _NDISC++;
       } while(true);
       //L1300: format(I1,9X,2F10.5);

*/



  //      int ID=1;
  //     PREPARE TABLES FOR RECOILING NUCLEUS.
  //      if(_IDIST!=0) MOMENT(0,0.,0.,0,_VZC,_NCASC,1,1);
  if(_IDIST!=0)
    for(I=1; I<=_NCASC; I++)
      MOMENT(1,0.,0.,0,_VZC,a[I], *htmlStream);
  //comment G1

  if(_INPUT!=2) {
      _MAXC= min(Max_MOM,_MAXC);

      for(I=2; I<=_MAXC; I++) _ALEV[I]+=_ALEV[I-1];
      SLEV=max(0.0001,_ALEV[_MAXC]);
      for(I=1; I<=_MAXC; I++) _ALEV[I]/=SLEV;

      //     RUNNING SUM FOR J RANDOM SELECTION FORMED
      for(L=1; L<=_NCASC; L++) {        //      L  IS NO. OF CASCADE BEING CALCULATED

          for(I=1; I<_MAXC; I++)  if(RANF()<=_ALEV[I]) break;
          a[L]->init(_IZ,_IN,_ENERGY-ECLOSS*RANF(),I,(I-1)*_MDIR);
        }
    }
  else
    for(L=1; L<=_NCASC; L++)
      a[L]->init(_IZ,_IN,_ENERGY,_MAXC,(_MAXC-1)*_MDIR);
  //----------------------------------------

  Emin_wait=1e10;

  //     COMPOUND NUCLEUS PREPARATION FINISHED
  //int LK=0;
  EU=0;   SC=0;

  IA=LA=MA=MERG=LERG=II=0;
  // = new fusion_event();
  //  c = new fusion_event;

  c.init(0,0);
  csi.init(0,0);


  TRACK(0,0,0,0,0,0,0,0,0,*htmlStream);   // put all to 0

L100:
  TRACK(0,0,0,EU,EU,_ENERGY,1,SPROB,ECLOSS, *htmlStream);

  IAC=0;       c.Z=0;

  for(I=1; I<=_NCASC; I++) {  //L109
      if(a[I]->Ex <= 0.)continue;
      IA=a[I]->Z+a[I]->N;

      if(IA < IAC)continue;
      else if(IA > IAC || (IA==IAC && c.Z==0)) {
          IAC=IA;                   //     LOWEST LIMIT FOR A
          LA=I;                     //     LOOP TO FIND HIGHEST A FOR PARTICULAR Z
          c.init(a[I]);
          c.Ex=0.;
          _ERGMIN=9999.;             //      ERGMIN IS MINIMUM ENERGY FOR GIVEN A
        }

      MA=I;                              //     HIGHEST LIMIT OF A

      csi.init(a[I]);

      if( !csi.testZN(c) )continue;

      if(a[I]->Ex > c.Ex) {      //     FIND HIGHEST EX FOR HIGHEST A
          c.Ex=a[I]->Ex;

          LERG=I;           //     FOR THE ABOVE LA , MA -  LERG , MERG ARE INDICES OF HIGHEST
        }                 //       ENERGY LIMITS

      if(a[I]->Ex == c.Ex)  MERG=I;

      _ERGMIN = min(a[I]->Ex,_ERGMIN);
    }

  //===============================================================

  if( IAC<6 )goto L201;
  //
  //     WHEN ALL ENERGIES ARE ZERO  -  THE WHOLE JOB IS FINISHED
  //
  if(c.Z==0)goto L100;
  SC=0.5*(c.Ai()%2);
  //     SC IS SPIN INDEX OF NUCLEUS -  0 OR 1/2

  if(Emin_wait  > c.Ex )
    {
      if(progress_inc == 0)
        {
          countMax = c.Ex;
          progress.setMaximum(countMax);
          progress_inc++;
          progress.setLabelText(QString::number(Emin_wait));
        }
      if(progress.wasCanceled()){
          return 10;
        }
      if(countMax-c.Ex > val)
        {
          val = countMax-c.Ex;
          progress.setValue(val);
          progress.setLabelText("<h3 style=\"color:#e80a20;\"> &#x420; &#x430; &#x0441; &#x0447; &#x0451; &#x0442; ... </h3>"
                                "<h3 style=\"color:#0f67e0;\">Calcul ...</h3><h3 style=\"color:#d0b534;\">"
                                "Calculando ...</h3><h3> current Excitation Energy = " + QString::number(c.Ex) + " MeV</h3>");
          qApp->processEvents();
        }
      // //FormWait->InputEx(Emin_wait);
    }

  for(I=1; I<=4; I++) { //L367
      II=5-I;
      AMASS(_RMASS,c,II,_MEBIN[II],_MAXJ[II],_MAXJS[II] ,_EBIN[II],_RLEV[II-1],IXPR[II],_BE, *htmlStream);
    }

  _VQS[1]=0.;
  if( (t=_IZR[2]+_INR[2]) > 0) _VQS[2]=1.44*_IZR[2]/(_RCCC[2]*pow(t,.3333));
  else                         _VQS[2]=0;

  if( (t=_IZR[3]+_INR[3]) > 0) _VQS[3]=2.88*_IZR[3]/(_RCCC[3]*pow(t,.3333));
  else                         _VQS[3]=0;




  //     GET MASS EXCESSES OF NUCLEI
  //
  //   RMASS  ARE MASSES OF RESIDUAL NUCLIDES
  //  THE EMIN ARE THE SMALLEST ENERGIES FOR WHICH RELIABLE CROSS SECTIONS
  //          ARE CALCULABLE.  THE FORMULAE WERE OBTAINED FROM EXPERIENCE.
  //
  //     EX OF RSID.NUC. - MAXIMUM   MAXIMUM EXCITATION ENERGY OF
  //     PARTICULAR RESIDUAL NUCLEUS

  for(I=1; I<=NODES; I++)   EXRSM[I]=c.Ex+_BE[I]-_EMIN[I];


  EXRSM[4]=c.Ex;



  for(K=1; K<=MODES; K++)
    {         //     FOR  EX.LT.0.  MODE NOT AVAILABLE

      if(EXRSM[K]<=0. || EXRSM[K]<=_EBIN[K][_MEBIN[K]] ) continue;
      if(_IZR[K]*(_IZR[K]*_INR[K]) >0) {
          fprintf(f09,"\\par Level density able for Z=%d  N=%d  does not go tohigh enough energy %.3f reduired"
                      "\\par MODE=%d  EBIN = %.3f  MEBIN = %d",
                  _IZR[K],_INR[K],EXRSM[K], K, _EBIN[K][_MEBIN[K]],_MEBIN[K]);
          *htmlStream << "<p>Level density able for Z = " << QString::number(_IZR[K]) << "  N = " << QString::number(_INR[K])
                      << " does not go to high enough energy " << QString::number(EXRSM[K]) << "reduired </p>"
                      << "<p> MODE = " << QString::number(K) << " EBIN = " << QString::number(_EBIN[K][_MEBIN[K]]) << "MEBIN = " << QString::number(_MEBIN[K])
                                                                                                                   << "</p>";
        }
    }
  htmlStream->flush();
L11:

  c.J=a[LERG]->J;           //     c.J  IS CASCADE SPIN INDEX


  for( I=1; I<=NODES; I++)
    if(EXRSM[I]>0.)goto L129;


  //
  //      EVAPORATION OF PARTICLES OF PARTICULAR Z,N WITH ENERGY E OR LESS NO
  //      LONGER POSSIBLE
  //
  if(EXRSM[4]>0.) IBOUND=1;
  else            goto L100;                  //  NEW Z,N OF MAX A NEEDED


L129:



  for(K=1; K<=MODES; K++) {
      // A AND EX IN PREDETERMINED LIMITS. J FIXED. PARTICLES
      // AND GAMMAS CAN BE EMITTED.
      if(EXRSM[K]<=0.) {
          IXRSM[K]=0;
          continue;
        }

      FIND(EXRSM[K],IXRSM[K],_EBIN[K], _MEBIN[K],K, *htmlStream);

      //     THE FOLLOWING IS A PROCEEDURE TO FIND THE INTEGER INDEX
      //     ASSOCIATED WITH A GIVEN EXCITATION ENERGY IN THE LEVEL
      //     FILE. EBIN IS THE ENERGY OF THE BIN IXRSM.

      EXMIN=EXRSM[K]-_EMAX[K]+_EMIN[K]-1.;
      FIND(EXMIN,IXMIN[K],_EBIN[K], _MEBIN[K],K, *htmlStream);
    }


L29:       //  CIKL here

  // FormWait->InputJ(c.J);

  //     ZERO THE PROBABILITIES BEFORE CALCULATION
  for(J=1; J<=5555; J++) _PROB[J]=0.;


  Fprob=0.;
  SPROB=0.;
  _IPROB=0;
  //     IPROBP COUNTS PARTICLE OUTPUT
  for(K=1; K<=MODES; K++)
    {
      PROBM[K]=0.;                      //     PROBM IS MAXIMUM PROBABILITY FOR MODE (CHANNEL) K.
      if(EXRSM[K]<=0.) continue;
      //     CHNPRB GIVES PROBABILITY FOR ALL CELLS IN MODE K
      CHNPRB(K,IXRSM[K],c,SC,_SPART[K],
             _RLEV[K-1],PROBM[K],_MAXJ[K],_MAXJS[K],_EBIN[K],IXMIN[K],*htmlStream);
    }


  if(_IPROB!=0)
    {
      //    AN EXCITED TRAP
      Fprob=_PROB[_IPROB]*_LFISS;
      for(I=2; I<=_IPROB; I++)    _PROB[I]+=_PROB[I-1];

      SPROB=_PROB[_IPROB];
      if(SPROB==0.)  _IPROB=0;
      else      {
          TPROB=1./SPROB;
          for(I=1; I<=_IPROB; I++) _PROB[I]*=TPROB;
          Fprob*=TPROB;
        }
    }

  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw start cikl L105
  for(I=LERG; I<=MERG; I++)
    {           //L105
      //     SCAN ALL EVENTS LOOKING FOR ONES WITH RIGHT A,Z,EX,J
      csi.init(a[I]);

      if( !csi.testZN(c)  || a[I]->Ex!= c.Ex ) continue;

      csi.J  = a[I]->J;
      csi.MJ = a[I]->MJ;

      if(csi.J != c.J) {
          if(LERG==0)LERG=I;
          continue;
        }         //     FOR THE FIRST DIFFERENT SPIN SET NEW LIMITS SO THAT IN THE NEXT
      //     SCNA WE USE SAME TL'S.

      if(LERG==I)LERG=0;

      //==========================================
      if(_IPROB==0)
        {           //     FOLLOWING INSERT FORCES DECAY THRU E2 GAMMA YRAST CASCADE

          ERGNOT=a[I]->Ex;
L3882:
          if(a[I]->J<=1)
            {
              a[I]->Ex=0;
              a[I]->J=0;
              csi.J=0;
              continue;
            }

          AJN=a[I]->J - 1 + SC;
          MODE=4;
          ETRAN=a[I]->Ex*(4.*AJN-2.)/(AJN*(AJN+1.));
          c.Ex=a[I]->Ex;
          a[I]->Ex=max(a[I]->Ex-ETRAN,0.);

          if(a[I]->J==c.J) {
              INOTGS++;
              ERGNOT=c.Ex;
            }
          csi.J=a[I]->J;
          a[I]->J-=2;
          //     NOT DECAYED TO G.S. (''YRAST TRAPS''- NO E3 OR M3 IN PROGRAM)
          //     END YRAST CASCADE INSERT.  15 MARCH 1979
          //==========================================
        }
      else  {
          //==========================================

          SEARCH(a[I],c,MODE, *htmlStream);

          if(a[I]->Ex<0) {    // fission
              csi.Ex = a[I]->Ex = 0.;
              TRACK(5,csi.J,csi.J,c.Ex,csi.Ex  ,c.Ex,2,SPROB,ECLOSS, *htmlStream);

              if(_IDIST!=0)
                {
                  //qDebug() << "HIGH_WORD(I) = " << HIGH_WORD(I);
                  IBUF2[0]=HIGH_WORD(I);
                  IBUF2[1]=LOW_WORD(I);
                  IBUF2[2]=(short int)6;
                  //     MODE=6 IS FISSION TRACEBACK
                  IBUF2[3]=(short int)csi.J;
                  IBUF2[4]=(short int)0;
                  IBUF2[5]=(short int)csi.MJ;
                  IBUF2[6]=(short int)(c.Ex+1.);
                  IBUF2[7]=(short int)0;
                  IBUF2[8]=(short int)c.Z;
                  IBUF2[9]=(short int)c.N;
                  //MK                    buf2List << HIGH_WORD(I) << LOW_WORD(I) << 6 << csi.J << csi.MJ << (c.Ex+1.) << 0 << c.Z << c.N;
                  //MK                    buf2Link << buf2List;
                  //MK                    buf2List.clear();
                  fwrite(IBUF2,sizeof(short int),10,f02);

                  BUF3[0]=(float)Fprob;
                  BUF3[1]=(float)c.Ex;
                  BUF3[2]=(float)a[I]->Ex;
                  BUF3[3]=(float)a[I]->particle_energy;
                  BUF3[4]=(float)a[I]->particle_angle;
                  //MK                    buf3List << Fprob << c.Ex << a[I]->Ex << a[I]->particle_energy << a[I]->particle_angle;
                  //MK                    buf3Link << buf3List;
                  //MK                    buf3List.clear();
                  fwrite(BUF3,sizeof(float),5,f02);
                }

              a[I]->Z=-a[I]->Z;
              //     NEGATIVE IZCS(I) IMPLIES FISSION OF IZCS(I)
              continue;
              //     SEARCH DETERMINES FINAL NUCLEUS
              //     PROCEED TO CATALOGUE DISTRIBUTIONS.
            }
        }
      //==========================================

      EP=c.Ex-a[I]->Ex+_BE[MODE];
      if(_IPROB==0)EP=ETRAN;
      csf.J=a[I]->J;
      csi.Ex=a[I]->Ex;
      TRACK(MODE,csi.J,csf.J, c.Ex,csi.Ex,c.Ex,2,SPROB,ECLOSS, *htmlStream);
      AC=csi.Ad();
      AP=_IZPART[MODE] + _INPART[MODE];
      EP*=AC/(AP+AC);
      int IE12=0;

      if(MODE==4) {
          int JDDJ=abs(csf.J-csi.J);

          if(JDDJ==2) IE12=1;
          else {
              double GMSUM=_GE1+_GE2*EP*EP;
              GMSUM=_GE1/GMSUM;
              //     IF DELTA-J IS 2 DECAY IS E2
              //     OTHERWISE, DECAY MODE (FOR TRACEBACK PURPOSES ONLY)
              //     IS DECIDED BY SIZE OF DECAY PROBABILITY (E1 OR E2 ONLY)
              double X=RANF();
              if(X>GMSUM)IE12=1;
            }
        }
      //      IBUF2(2) NOW GIVES MODE=4 FOR E1 GAMMAS AND MODE=5 FOR E2 GAMMAS
      //      *** THIS IS FOR THE DUMP-TRACEBACK ONLY

      if(_IDIST!=0)
        {
          // qDebug() << "HIGH_WORD(" << I << ") = " << HIGH_WORD(I);
          //qDebug() << "LOW_WORD(" << I << ") = " << LOW_WORD(I);
          IBUF2[0]=HIGH_WORD(I);
          IBUF2[1]=LOW_WORD(I);
          IBUF2[2]=(short int)(MODE+IE12);

          if(MODE+IE12 > 4) EP+=0.00001;

          IBUF2[3]=(short int)csi.J;
          IBUF2[4]=(short int)csf.J;
          IBUF2[5]=(short int)csi.MJ;
          IBUF2[6]=(short int)(c.Ex+1.);
          IBUF2[7]=(short int)(EP*10.+1.);
          IBUF2[8]=(short int)c.Z;
          IBUF2[9]=(short int)c.N;
          //MK            buf2List << csi.J << csf.J << csi.MJ << (c.Ex + 1.) << (EP+10.+1.) << c.Z << c.N;
          fwrite(IBUF2, sizeof(short int), 10, f02);
          //MK            buf2Link << buf2List;
          //MK            buf2List.clear();
        }


      int JKQ=EP+1.;

      JKQ = min(JKQ,50);

      NSPC[MODE][JKQ]++;

      if(_IDIST!=0)
        {
          if(MODE==MODES)
            {
              int JCSFF=csf.J-1;
              a[I]->MJ = min(max(a[I]->MJ,-JCSFF), JCSFF);
              a[I]->particle_energy = EP;
              OUTEM((IBOUND==1 ? 4 : 2),*htmlStream,MODE,EP);                   //     STORE GAMMA ENERGY
              a[I]->particle_angle = 0;
            }
          else    {
              csf.MJ=0;
              MJRAN(csi,csf,MODE,EP, *htmlStream);
              a[I]->MJ=csf.MJ;
              MOMENT(2,csi.Ad(),AP,MODE,EP,a[I],*htmlStream);               //     STORE PARTICLE ENERGY. FORMS OUTPUT DISTRIBUTIONS BYING OUTEM().
            }
        }

      if(_IDIST!=0)
        {
          BUF3[0]=(float)Fprob;
          BUF3[1]=(float)c.Ex;
          BUF3[2]=(float)a[I]->Ex;
          BUF3[3]=(float)a[I]->particle_energy;
          BUF3[4]=(float)a[I]->particle_angle;

          fwrite(BUF3,sizeof(float),5,f02);
        }
      if(a[I]->Ex>0.&&_IPROB==0) goto L3882;
    }

  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw end cikl L105



  if(_IPROB==0) {
      ENOTGS += ERGNOT*(INOTGS-INOTG1);     // for Average energy
      SNOTGS += c.J*   (INOTGS-INOTG1);     // for average spin
      INOTG1=INOTGS;
    }

  if(LERG!=0) {               //     NO FURTHER SPIN AT THIS A AND Z
      c.J = a[LERG]->J;   //     NEW ENERGY OR NUCLIDE NEEDED
      goto L29;          //     NEW SPIN NEEDED
    }


  IBOUND=0;
  c.Ex=0.;

  for( I=LA; I<=MA; I++) {
      if( a[I]->Ex == 0     ||
          a[I]->Ex < c.Ex   ||
          !c.testZN(a[I])  )   continue;

      if(a[I]->Ex>c.Ex) {
          c.Ex=a[I]->Ex;
          LERG=I;
        }
      MERG=I;
    }



  if(c.Ex>0) {
      for(K=1; K<=NODES; K++) EXRSM[K] = c.Ex + _BE[K] - _EMIN[K];
      EXRSM[4]=c.Ex;
      goto L11;
    }

  goto L100;
  //    NEW NUCLIDE OF SAME A NEEDED

  //===================================  FINISH
L201:
  htmlStream->flush();
  if(_IDIST!=0)
    for(I=1; I<=_NCASC; I++)
      if(a[I]->Z>0) MOMENT(3,a[I]->Ad(),0.,0,0.,a[I], *htmlStream);

  //     LAST PARAMETER IN MOMENT
  //        = 1  INITIALIZATION OF TABLES
  //        = 2  COLLECTING DATA DURING RUN
  //        = 3  TERMINATION AND OUTPUT.


  fprintf(f09,"\\par\\par\\par\\f2\\cf3\\b\\fs30   ------   Output results for compound nucleus decay -----\\par\\f0\\b0\\cf0\\fs16");
  *htmlStream << "<p>&nbsp;</p><p>&nbsp;</p><h2 align=center> Output results for compound nucleus decay </h2>";
  II=1;
  htmlStream->flush();
  PRODCT(II,IZMAX,INMAX,_NCASC,a,filename_cs, *htmlStream);
  if(_IDIST!=0) OUTEM(3,*htmlStream);
  if(INOTGS!=0) {
      ENOTGS/=double(INOTGS);
      double PR=(100.*INOTGS)/_NCASC;
      SNOTGS/=double(INOTGS);
      fprintf(f09,
              "\\fs20\\f2\\i\\fc2\\cf2 "
              "\\par\\par\\par  \\i0\\b %.2f\\i\\b0  percent of cascades trapped before reaching ground state due to spin inhibition"
              "\\par  Average energy at which cascades were trapped is \\i0\\b %.2f\\i\\b0  MeV, average spin =\\i0\\b %.1f\\i\\b0  hbar"
              "\\par  **** successive decays through single yrast cascade assumed"
              "\\fs16\\f0\\i0\\cf0"
              ,PR,ENOTGS,SNOTGS);
      *htmlStream << "<p>&nbsp;</p><p>&nbsp;</p><p><b>" << QString::number(PR) << "  percent of cascades trapped before reaching ground state due to spin inhibition</b></p>"
                                                                                  "<p> Average energy at which cascades were trapped is <b><em>" << QString::number(ENOTGS) << "</em></b>   MeV, average spin = <b><em>"
                  << QString::number(SNOTGS) << "</em></b> hbar </p><p> **** successive decays through single yrast cascade assumed </p>";
    }

  fprintf(f09,"\\par"
              "\\par\\b --------  C.M. spectra  of  emitted  particles ------\\b0 "
              "\\par\\par\\b\\cf1        Ex(MeV)     Neut     Prot    Alpha    Gamma   \\b0\\cf0 ");
  *htmlStream << "<p>&nbsp;</p><h3> --------  C.M. spectra  of  emitted  particles ------ </h3><p>&nbsp;</p>"
                 "<table border=1 cellspacing=0 cellpadding=2><tr><th>Ex(MeV)</th><th>Neut</th><th>Prot</th><th>Alpha</th><th>Gamma</th></tr>";


  for(IST=1; IST<=4; IST++) _BE[IST]=0.;

  for(int IKQ=1; IKQ<=50; IKQ++) {
      int IKQ1=IKQ-1;
      int IKZ=IKQ;
      if(IKQ==50)IKZ=99;
      if(NSPC[1][IKQ]!=0 ||
         NSPC[2][IKQ]!=0 ||
         NSPC[4][IKQ]!=0 ||
         NSPC[3][IKQ]!=0)
        {
          fprintf(f09,"\\par |    %3d - %-3d  |",IKQ1,IKZ);
          *htmlStream << "<tr><td align=center>" << QString::number(IKQ1) << " - " << QString::number(IKZ) << "</td>";
          for(I=1;I<=4;I++)
            if(NSPC[I][IKQ]>0) {fprintf(f09," %6d |",NSPC[I][IKQ]);
                *htmlStream << "<td align=center>" << QString::number(NSPC[I][IKQ]) << "</td>";
              }else   {            fprintf(f09,"        |");
                *htmlStream << "<td></td>";
              }
          *htmlStream << "</tr>";
        }
      double EPART=double(IKQ)-0.5;
      for(IST=1; IST<=4; IST++) {
          NSPC[IST][51]+=NSPC[IST][IKQ];
          _BE[IST]+=EPART*NSPC[IST][IKQ];
        }
    }

  for(IST=1; IST<=4; IST++) _BE[IST]/=(NSPC[IST][51]+1.E-9);


  fprintf(f09,"\\par |   \\b   Total   \\b0   |\\b  %6d \\b0 |\\b  %6d \\b0 |\\b  %6d \\b0 |\\b  %6d \\b0 |",
          NSPC[1][51],NSPC[2][51],NSPC[3][51],NSPC[4][51]);
  *htmlStream << "<tr><th> Total  </th><td align=center>" << QString::number(NSPC[1][51]) << "</td><td align=center>" << QString::number(NSPC[2][51]) << "</td><td>"
                                                                                                                                                      << QString::number(NSPC[3][51]) << "</td><td> " << QString::number(NSPC[4][51]) << "</td></tr>";
  fprintf(f09,"\\par | AverageEnergy | %6.2f | %6.2f | %6.2f | %6.2f |",
          _BE[1], _BE[2], _BE[3], _BE[4]);
  *htmlStream << "<tr><th> Average Energy</th><td align=center>" << QString::number(_BE[1]) << "</td><td align=center> " << QString::number(_BE[2]) << "</td><td align=center>"
                                                                                                                                                    << QString::number(_BE[3]) << "</td><td align=center>" << QString::number(_BE[4]) << "</td></tr></table>";
  fprintf(f09,"\\par -----------------------------------------------------\\par ");
  //*htmlStream << "<p> -----------------------------------------------------</p>";
  TRACK(0,0,0,c.Ex,c.Ex,c.Ex,3,SPROB,ECLOSS, *htmlStream);

  fprintf(f09,"\\par\\par\\b\\cf2     ---- end of evaporation calculation ----\\b0\\cf0 ");
  *htmlStream << "<p>&nbsp;</p><p>     ---- end of evaporation calculation ---- </p>";

  if(_IDIST!=0 && f02)  {
      fclose(f02); f02=nullptr;
      STATIS(a, filename_evt, ParticleFlags, *htmlStream);
    }

  if(f_particles)  {
      fclose(f_particles);
      f_particles=nullptr;
      ParticleFlags.ClearParticleState(ep_opened);
    }
  *htmlStream << "</body></html>";
  htmlStream->flush();
  // progress.setValue(countMax);
  // qApp->processEvents();
BadEnd:
  fprintf(f09,"}");
  fclose(f09);
  htmlFile->close();

  for(I=0; I<=_NCASC; I++) delete a[I];
  delete a;

  return 1;
}




