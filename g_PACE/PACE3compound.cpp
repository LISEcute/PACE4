#include <QDebug>
#include "WinLise_Constant.h"
#include "ParticleStates.h"
#include "ftype.h"
#include <math.h>
#include <stdio.h>
#include "pace.h"
#include <QTextStream>

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))


void TLLL(double ER, int NP, int IZ, int IN, double TL[25], int &LLMAX, QTextStream &s);
extern double RANF();
extern double mzsqrt(double X);
extern double pow2(double par);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
extern FILE *f09;

extern const char *s40star;
extern const char *s10star;
//extern char *s35x;
const char s35x[10] = "         ";
extern char *s100m;


///-------------------- static to clear
extern int _cl_K;    ///   it is necessary to clear each time  CHECK before run PACE3 !!!!!!
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::COMPND(int IZ, int IN, double ENERGY, int MAXC, QTextStream& s)
{
  //
  //     ORIGINAL JULIAN SUBROUTINE
  //
  //     ***********
  //     FORM J DISTRIBUTION OF INITIAL COMPOUND NUCLEI

  //=========      common/CP/
  extern double _ALEV[Max_MOM+1];


  double AJ, D, AL,FL,SIG,F1,ACONST, UEX;

  double Z=IZ;
  double A=IZ+IN;
  double SR=((IZ+IN)%(2))*.5;

  if(A*Z > 0)
    {
      GCLDP(A,Z,0,FL,AL,SIG,F1,D,ACONST,SR,1.,s);


      UEX=ENERGY-D;
      SIG=1./(2.*SIG*sqrt(UEX/AL));

      if(MAXC>Max_MOM)MAXC=Max_MOM;

      for(int J=1; J<=MAXC; J++) {
          AJ=J-1+SR;
          _ALEV[J]=(2.*AJ+1.)*exp(-pow((AJ+.5),2)*SIG);
        }
    }
  else _ALEV[1]=1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::GCLDP (double AMR, double AZR, int IADEF, double &FALIT,
                  double &ALIT, double &SIGSQ, double &FACT1,
                  double &DELTA, double &ACONST, double &SR, double ARATIO, QTextStream& s)
{

  //==================     common /XQTT/
  extern double _UX,_EX,_TT;
  //=========      common/QMMJ/
  extern double _FACLA,_FYRST,_EROT[Max_MOM+1];
  extern int    _MAXC;


  //***  LEVEL DENSITY PARAMETERS WERE TAKEN FROM A. GILBERT AND A.G.W.
  //***  CAMERON, CANADIAN JOURNAL OF PHYSICS, VOLUME 43 (1965) 1446.
  //***  DATA DO NOT INCLUDE LEVEL DENSITY PARAMETERS FOR Z OR N LT.11

  // base 1 / already increased
  double PZGC[99] = {
    0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.46,0.,2.09,0.,1.62,0.,1.62,0.,1.83,0.,1.73,0.,1.35,
    0.,1.54,0.,1.20,0.,1.06,0.,1.36,0.,1.43,0.,1.17,0.,1.24,0.,1.20,0.,1.28,0.,1.28,0.,1.35,
    0.,1.36,0.,1.19,0.,1.14,0.,1.12,0.,1.58,0.,1.17,0.,1.18,0.,1.22,0.,0.97,0.,0.92,0.,0.62,
    0.,0.68,0.,0.64,0.,0.72,0.,0.75,0.,0.71,0.,0.87,0.,0.83,0.,0.89,0.,0.79,0.,0.89,0.,0.78,
    0.,0.69,0.,0.61,0.,0.72,0.,0.77};
  double PNGC[151] = {
    0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.67,0.,1.80,0.,1.67,0.,1.86,0.,2.04,0.,1.64,0.,1.44,
    0.,1.54,0.,1.30,0.,1.27,0.,1.29,0.,1.41,0.,1.50,0.,1.50,0.,1.43,0.,1.88,0.,1.47,0.,1.57,
    0.,1.46,0.,0.93,0.,0.72,0.,1.12,0.,1.29,0.,0.94,0.,1.24,0.,1.25,0.,1.14,0.,1.32,0.,1.15,
    0.,1.24,0.,1.43,0.,1.09,0.,1.20,0.,1.04,0.,0.70,0.,0.85,0.,0.76,0.,0.92,0.,0.99,0.,1.10,
    0.,0.92,0.,0.73,0.,0.70,0.,0.87,0.,0.61,0.,0.69,0.,0.55,0.,0.40,0.,0.73,0.,0.58,0.,0.86,
    0.,1.13,0.,0.84,0.,0.79,0.,0.82,0.,0.71,0.,0.41,0.,0.38,0.,0.67,0.,0.61,0.,0.78,0.,0.67,
    0.,0.67,0.,0.79,0.,0.60,0.,0.57,0.,0.49,0.,0.43,0.,0.50,0.,0.39};
  double SZGC[99] = {
    0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-2.91,-4.17,-5.72,-7.80,-8.97,-9.70,-10.10,-10.70,-11.38,
    -12.07,-12.55,-13.24,-13.93,-14.71,-15.53,-16.37,-17.36,-18.52,-18.44,-18.19,-17.68,-17.09,
    -16.65,-16.66,-16.59,-16.35,-16.18,-16.41,-16.60,-16.54,-16.42,-16.84,-17.22,-17.42,-17.52,
    -17.82,-18.19,-18.58,-19.11,-19.83,-19.14,-18.35,-17.40,-16.54,-15.68,-14.75,-13.71,-12.87,
    -12.18,-11.61,-11.09,-10.78,-10.53,-10.41,-10.21,-9.85,-9.36,-8.97,-8.56,-8.13,-7.68,-7.33,
    -7.11,-7.16,-7.05,-6.81,-6.56,-6.95,-7.52,-8.03,-8.41,-8.86,-7.71,-6.38,-5.47,-4.78,-4.37,
    -4.17,-4.12,-4.29,-4.61,-5.04,-5.48,-5.96,-6.40,-6.87,-7.20,-7.74};
  double SNGC[151] = {
    0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,6.80,7.53,7.55,7.21,7.44,8.07,8.94,9.81,10.60,11.39,12.54,
    13.68,14.34,14.19,13.83,13.50,13.00,12.13,12.60,13.26,14.13,14.92,15.60,16.38,17.08,17.55,
    17.98,18.33,18.56,18.71,18.65,18.55,18.52,18.34,18.01,17.38,16.56,15.62,14.38,12.88,13.24,
    13.71,14.40,15.16,15.89,16.43,16.97,17.59,18.08,18.72,19.22,19.51,19.73,19.91,20.06,20.16,
    20.09,19.83,19.41,19.06,18.66,17.73,17.03,16.44,16.00,15.33,14.49,13.42,12.28,11.14,10.10,
    9.09,10.00,10.64,11.18,11.70,12.22,12.71,13.05,12.99,12.62,12.11,11.66,11.21,10.81,10.38,
    10.03, 9.65, 9.38, 8.99, 8.62, 8.33, 8.10, 7.82, 7.56, 7.33, 7.15, 6.83, 6.69, 6.55, 6.53,
    6.49, 6.39, 5.82, 5.26, 4.53, 3.83, 3.08, 2.37, 1.72, 1.05, 0.27,-0.69,-1.69,-2.58,-3.16,
    -1.72,-0.41, 0.71, 1.66, 2.62, 3.22, 3.76, 4.10, 4.46, 4.83, 5.09, 5.18, 5.17, 5.10, 5.05,
    5.04, 5.03, 4.99, 4.98, 5.11, 5.27, 5.39, 5.37, 5.30};

  int IZR,INR,MAXXL;
  double CNR,AL;


  IZR=AZR+0.01;
  CNR=AMR-AZR+0.5;

  CNR=min(150,CNR);    // Oleg   -  delete after extrapolation
  INR=CNR;
  AZR=min(AZR,98);      // Oleg   -  delete after extrapolation

  if(_FACLA!=0.)       ALIT=AMR/_FACLA*ARATIO;
  else {                        //  GILBERT & CAMERON BYPASSED FOR LITTLE-A
      if( IZR<11 || INR<11 ) {
          fprintf(f09,"\\par Level density parameters for Z=%-3d and/or  N=%-3d are "
                      "not included as double in void GCLDP",IZR,INR);
          s << "<p> Level density parameters for Z = " << QString::number(IZR) << " and/or N = " << QString::number(INR)
            <<  "not included as double in void GCLDP</p>";
          ACONST=0.;
          s.flush();
          return;
        }

      ACONST=(IADEF==2 ? 0.142 : 0.120);
      FALIT=0.00917*(SZGC[IZR]+SNGC[INR])+ACONST;
      ALIT = FALIT * AMR * ARATIO;
    }

  // 05/27/2010 Oleg ALIT should be greater than 0

  SIGSQ=0.0888*ALIT*pow(AMR,0.66667);

  FACT1=12.0*pow(8.0,0.5)*pow(ALIT,0.25)*pow(SIGSQ,1.5);
  if(FACT1==0)
    FACT1=1;
  FACT1=1./FACT1;

  //      OBTAIN ROTATION ENERGIES FROM C.P.S.
  MAXXL=min(Max_MOM , _MAXC+9);

  for(int L=1; L<=MAXXL; L++) {
      AL=L-1+SR;
      if(_FYRST!=0.) { _EROT[L]=YRAST(AMR,AZR,AL)*_FYRST;  break;}
      _EROT[L]=(AL+0.5)*(AL+0.5)/(2.*SIGSQ);
    }


  DELTA=PZGC[IZR]+PNGC[INR];
  _UX=2.5+150./AMR;
  _EX=_UX+PZGC[IZR]+PNGC[INR];
  _TT=sqrt(ALIT/_UX)-1.5/_UX;
  _TT=1./_TT;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double PACE::YRAST(double A, double Z, double AL) {

  //     ROTATING LIQUID DROP YRAST LINE.
  //     ***** COURTESY OF DR. F. PLASIL, O.R.N.L.

#include "yrast.h"   // inverse array!!!

  double AN=A-Z;
  double PAREN=1.-1.7826*pow(((AN-Z)/A),2);
  double ESO=17.9439*PAREN*pow(A,0.66667);
  double X=0.019655*Z*(Z/A)/PAREN;
  double ERO=34.548*(AL*AL)/pow(A,1.6667);
  double Y=1.9254*(AL*AL) /(PAREN*pow(A,2.3333));
  int    IX=20.*X+1.;
  double CX=IX;
  double BX=20.*X+1.;
  double DX=BX-CX;
  double BY, CY,DY, H1, H2, HF;
  int IY;

  if(X<=0.25)
    {
      BY=10.*Y+1.;
      BY=min(BY,9);    BY=max(BY,1);
      IY=BY;  CY=IY;   DY=BY-CY;
      H1=(X1H[IY  ][IX+1]-X1H[IY  ][IX])*DX+X1H[IY  ][IX];
      H2=(X1H[IY+1][IX+1]-X1H[IY+1][IX])*DX+X1H[IY+1][IX];
    }
  else
    {
      if(X<=0.5)
        {
          BY=20.*Y+1.;
          BY=min(BY,11);   BY=max(BY,1);
          IX-=5;
          IY=BY; CY=IY;  DY=BY-CY;
          H1=(X2H[IY  ][IX+1]-X2H[IY  ][IX])*DX+X2H[IY  ][IX];
          H2=(X2H[IY+1][IX+1]-X2H[IY+1][IX])*DX+X2H[IY+1][IX];
        }
      else
        {
          if(X>0.95)X=0.95;
          IX=20.*X+1.;
          IX-=10;
          BY=100.*Y+1.;
          BY=min(BY,19);  BY=max(BY,1);
          IY=BY;  CY=IY;  DY=BY-CY;
          H1=(X3H[IY  ][IX+1]-X3H[IY  ][IX])*DX+X3H[IY  ][IX];
          H2=(X3H[IY+1][IX+1]-X3H[IY+1][IX])*DX+X3H[IY+1][IX];
        }
    }

  HF=(H2-H1)*DY+H1;
  return ERO+HF*ESO;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double PACE::FISROT(double A, double Z, double AL, double BARFAC) {
  //     COHEN-PLASIL-SWIATECKI ROTATING LIQUID DROP FISSION BARRIER
  //     PROVIDED BY COURTESY OF DR. F. PLASIL, OAK RIDGE NATIONAL LAB.
  //
  //     BARRIER=DELSP-DELR
  //
  // it was done like here
  //  http://wwwasd.web.cern.ch/wwwasd/geant/geant4_public/G4UsersDocuments/
  //  UsersGuides/PhysicsReferenceManual/html/node62.html#SECTION031142000000000000000

#include "fisrot.h"         // inverse array!!!   X1B , X2B, X3B   [X] <-> {Y]

  double BF;
  double BY,  DY, B1, B2;
  int IY;


  const double  k = 1.7826;
  const double aS = 17.9439; // MeV
  const double aC =  0.7053; // MeV

  double N=A-Z;
  double PAREN=1.-k*pow2((N-Z)/A);        // k =1.7826  Cohen S. and Swiatecki W. J., Ann. Phys. 22 406 (1963).
  double X=(aC/2/aS)*(Z*Z/A)/PAREN;       // the fissility parameter

  double ESO=aS*PAREN*pow(A,2./3.);
  double Y=1.9254*AL*AL/(PAREN*pow(A,2.3333));

  double BX=20.*X+0.999;
  int    IX=BX;
  double DX=BX-IX;

  if(X<=0.25) {
      BY=10.*Y+0.999;
      BY = min(BY,9);  BY=max(BY,1);
      IY=BY;  DY=BY-IY;
      B2=(X1B[IY+1][IX+1] - X1B[IY+1][IX])*DX + X1B[IY+1][IX];
      B1=(X1B[IY  ][IX+1] - X1B[IY  ][IX])*DX + X1B[IY  ][IX];
      BF=(B2-B1)*DY+B1;
    }
  else       {
      if(X<=0.5) {
          BY=20.*Y+0.999;
          BY = min(BY,11);  BY=max(BY,1);
          IX-=5;
          IY=BY;  DY=BY-IY;
          B2=(X2B[IY+1][IX+1] - X2B[IY+1][IX])*DX + X2B[IY+1][IX];
          B1=(X2B[IY  ][IX+1] - X2B[IY  ][IX])*DX + X2B[IY  ][IX];
          BF=(B2-B1)*DY+B1;
        }
      else       {
          X = min(X,0.95);
          IX=20.*X+0.999;
          IX-=10;
          BY=100.*Y+0.999;
          BY = min(BY,19);  BY=max(BY,1);
          IY=BY;  DY=BY-IY;
          B2=(X3B[IY+1][IX+1] - X3B[IY+1][IX])*DX + X3B[IY+1][IX];
          B1=(X3B[IY  ][IX+1] - X3B[IY  ][IX])*DX + X3B[IY  ][IX];
          BF=(B2-B1)*DY+B1;
        }
    }

  BF*=BARFAC*ESO;

  return BF;
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double PACE::SIERK(double A , double Z, double AL, double BARFAC) {

  double BRR, BRL, DL, F;
  int J1, J2, L1, K, L;
  // ARRAY has been inverted
  double BS[3][17] = {      // BASe = 0 !!!!!!!!!!!!!   Oleg
                            {25.15, 25.01, 24.63, 23.99, 23.11, 21.98, 20.63, 19.07, 17.29, 15.32, 13.19, 10.90, 8.49,  5.98,  3.42, 0.91,  0.00},
                            {18.61, 18.51, 18.22, 17.75, 17.08, 16.24, 15.24, 14.06, 12.74, 11.29,  9.72,  8.05, 6.30,  4.52,  2.77, 1.15,  0.00},
                            {10.51, 10.44, 10.24,  9.90,  9.44,  8.86,  8.17,  7.37,  6.50,  5.55,  4.56,  3.54, 2.53,  1.58,  0.79, 0.32,  0.00}
                     };


  //=========      common/QSIER/
  extern double _BARLD[81];
  extern int  _IS;
  //========================


  if(_IS!=0)
    {
      L=AL;

      if(L>=80) return 0;

      J1=(L+5)/5;
      J2=J1+1;
      L1=(J1-1)*5;
      DL=L-L1;

      J1--; J2--; K=_IS-1;  // Oleg gor BASE=0

      BRR=BS[K][J1]+(BS[K][J2]-BS[K][J1])*0.2*DL;

      BRL = FISROT(A,Z,AL,1.);
      F=BRL/(_BARLD[L+1]+1.e-9);
      BRR *= F * BARFAC;
    }
  else
    BRR = BARFAC*FISROT(A,Z,AL,1.);


  return BRR;
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::OUTEM(int ICTL, QTextStream& s, int MODE, double EEML, double AEML) {
  //
  //     ORIGINAL JULIAN SUBROUTINE
  //
  //     ***********
  //
  //   ICTL IS ENTRY CONTROL
  //   MODE IS KIND OF PARTICLE EMITTED OR GAMMA
  //   EEML IS LABORATORY ENERGY OF EMITTED PARTICLE OR GAMMA
  //   AEML IS LABORATORY ANGLE OF EMITTED PARTICLE
  //   NT IS NUMBER OF EVENTS OF EACH MODE
  //   N(MODE,ENERGY,ANGLE) IS DISTRIBUTION OF EVENTS
  //       31ST ENERGY IS FOR ENERGY OVERFLOW
  //       32ND ENERGY IS FOR TOTAL AT EACH ANGLE
  //       THERE ARE 18 ANGLES
  //   NG(ENERGY) IS FOR EACH ENERGY NUMBER OF GAMMAS
  //   NW(MODE,BIN,ANGLE) IS DISTRIBUTION OF EVENTS INTO FOUR ENERGY BINS
  //   NGW(BIN) IS DISTRIBUTION OF GAMMAS INTO FOUR ENERGY BINS
  //   EW IS ENERGY OF SEPARATION BETWEEN FOUR ENERGY BINS
  //   EW1 IS LOWER ENERGY OF EACH OF FOUR ENERGY BINS
  //   EW2 IS UPPER ENERGY OF EACH OF THREE ENERGY BINS
  //                       THE FOURTH IS INFINITY
  //


  //--------common /XQSIG/
  extern double _SIGMA;
  extern int _NCASC;
  extern double _DSIG[19];
  //==================     common/EM/
  extern int eNT[5],eN[4][33][19],eNG[33],eNW[4][5][19],eNGW[5],eNGB[33],eNGBW[5],eNGT[33],eNGTW[33];
  extern double eEW[4],eEW1[5],eEW2[4];


  //--------------
  double IWE=1;
  double DELE=1;
  int I,K,L,KW, KGM, ISUM, IQ;
  double AK, E1, E2, TET;
  char s500[130]= "\\par=============================================================================================================================";
  //char *s500html= "<p>=============================================================================================================================</p>";
  char s501[130]= "\\par|--------------|------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|";
  //char *s501html= "<p>|--------------|------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|</p>";
  char s500a[92]="\\par=======================================================================================";
  char s500ahtml[95]="<p>=======================================================================================</p>";
  char s502[92]= "\\par|------------------------------------|------------------|---------------|-------------|";
  //char *s502html= "<p>|------------------------------------|------------------|---------------|-------------|</p>";
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw

  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw

  switch(ICTL)   {

    //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw

    case 1 :
      //   DELE IS THE ENERGY DIFFERENCE FOR THE 30 ENERGY BINS
      //   IWE IS CONTROL  0 IS NO FOUR ENERGY BINS  1 IS YES
      //
      //   INITIATION
      //
      //     **************************

      eEW[1]=5.;
      eEW[2]=10.;
      eEW[3]=20.;

      for(I=1; I<=4; I++)  eNT[I]=0;

      for(K=1; K<=32; K++) {
          eNG[K]=eNGB[K]=eNGT[K]=0;
          for(I=1; I<=3; I++)
            for(L=1; L<=18; L++)
              eN[I][K][L]=0;
        }

      if(IWE==0)break;   // ??

      for( KW=1; KW<=4; KW++) {
          eNGW[KW]=eNGBW[KW]=eNGTW[KW]=0;
          for(I=1; I<=3; I++)
            for(L=1; L<=18; L++)
              eNW[I][KW][L]=0;
        }



      eEW1[1]=0.;
      eEW2[1]=eEW1[2]=eEW[1];
      eEW2[2]=eEW1[3]=eEW[2];
      eEW2[3]=eEW1[4]=eEW[3];
      break;
      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw

      //
      //   MONITORING
      //
    case 2:
    case 4:
      eNT[MODE]++;
      KGM =  min(EEML+1.0001,31);
      K =    min(EEML/DELE+1,31);
      //      ENERGY BINS

      if(MODE!=4) {
          L=AEML/10.+1;
          //      ANGLE BINS
          eN[MODE][K][L]++;
          eN[MODE][32][L]++;
        }
      else{
          eNGT[KGM]++;
          eNGT[32]++;
          if(ICTL!=4) {
              eNG[KGM]++;
              eNG[32]++;
            }
          else {
              eNGB[KGM]++;
              eNGB[32]++;
            }
        }

      if(IWE==0 || MODE==4) break;

      if(EEML<eEW[1])     KW=1;
      else if(EEML<eEW[2])     KW=2;
      else if(EEML<eEW[3])     KW=3;
      else                     KW=4;

      eNW[MODE][KW][L]++;
      break;
      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
      //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
      //
      //   GRAPHING
      //
    case  3:

      //----------------------------------------  start cikl L36
      for(MODE=1; MODE<=3; MODE++) {
          fprintf(f09,"\\par\\par\\par  %s\\b\\cf1\\f2 ",s35x);
          s << "<p>&nbsp;</p><p>";// << QString::fromUtf8(s35x) << "</p>";
          switch(MODE) {
            case 1 : fprintf(f09,"Neutron\\cf3  spectra in laboratory coordinates\\cf0  %d\\cf3  events\\cf0 ",eNT[1]);
              s << "<span style=\"color:red\">Neutron </span> spectra in laboratory coordinates " << QString::number(eNT[1]) << " events</p>";
              break;
            case 2 : if(eNT[2]==0)continue;
              fprintf(f09,"Proton\\cf3  spectra in laboratory coordinates\\cf0  %d\\cf3  events\\cf0 ",eNT[2]);
              s << "<span style=\"color:red\">Proton </span> spectra in laboratory coordinates " << QString::number(eNT[2]) << " events</p>";
              break;
            case 3 : if(eNT[3]==0)continue;
              fprintf(f09,"Alpha\\cf3  spectra in laboratory coordinates\\cf0  %d\\cf3  events\\cf0 ",eNT[3]);
              s << "<span style=\"color:red\">Alpha </span> spectra in laboratory coordinates " << QString::number(eNT[3]) << " events</p>";
              break;
            };

          fprintf(f09,"\\b0\\f0\\fs14 %s"
                      "\\par | Energy range |_________________________________________ Angular range (deg)_______________________________________________|"
                      "\\par |    (MeV)     |    0 |  10 |  20 |  30 |  40 |  50 |  60 |  70 |  80 |  90 | 100 | 110 | 120 | 130 | 140 | 150 | 160 | 170 |"
                      "\\par |              |   10 |  20 |  30 |  40 |  50 |  60 |  70 |  80 |  90 | 100 | 110 | 120 | 130 | 140 | 150 | 160 | 170 | 180 |"
                      "%s",s500,s501);

          s << "<table border=1 cellspacing=0 cellpadding=1><tr><th>Energy range</th><th colspan=18>Angular range (deg)</th></tr><tr><th>(MeV)</th><th>0</th><th>10</th><th>20</th><th>30</th><th>40</th><th>50</th><th>60"
               "</th><th>70</th><th>80</th><th>90</th><th>100</th><th>110</th><th>120</th><th>130</th><th>140</th><th>150</th><th>160</th><th>170</th></tr><tr><th></th>"
               "<th>10</th><th>20</th><th>30</th><th>40</th><th>50</th><th>60</th><th>70</th><th>80</th><th>90</th><th>100</th><th>110</th><th>120</th><th>130</th><th>140"
               "</th><th>150</th><th>160</th><th>170</th><th>180</th></tr>";
          for( K=1; K<=30; K++) {
              AK=K;
              E1=(AK-1.)*DELE;
              E2=E1+DELE;
              ISUM=0;
              for(IQ=1; IQ<=18; IQ++) ISUM+=eN[MODE][K][IQ];

              if(ISUM <= 0) continue;

              fprintf(f09,"\\par |%6.1f-%-6.1f | ",E1,E2);
              s << "<tr><td align=center>" << QString::number(E1) << "-" << QString::number(E2) << "</td>";
              //s.setFieldWidth(6);
              for(IQ=1; IQ<=18; IQ++){
                  if(eN[MODE][K][IQ]>0) {
                      fprintf(f09,"%5d|",eN[MODE][K][IQ]);
                      s << "<td>" << qSetFieldWidth(4) << QString::number(eN[MODE][K][IQ]) << "</td>";

                    }else{
                      fprintf(f09,"     |");
                      s << "<td>" << qSetFieldWidth(4) << "</td>";
                    }
                }
              s << "</tr>";
            }


          E1=30.*DELE;
          fprintf(f09,"\\par | Above %6.1f | ",E1);
          s << "<tr><td> Above " << QString::number(E1) << " </td>";
          for(IQ=1; IQ<=18; IQ++)
            if(eN[MODE][31][IQ]>0) {
                s << "<td>" << QString::number(eN[MODE][31][IQ]) << "</td>";
                fprintf(f09,"%5d|",eN[MODE][31][IQ]);
              }else  {                  fprintf(f09,"     |");
                s << "<td></td>";
              }
          s << "</tr>";
          fprintf(f09,"%s",s501);
          // s << QString::fromUtf8(s501html);
          E1=0.;

          fprintf(f09,"\\par | Total        | ");
          s << "<tr><td>Total</td>";
          for(IQ=1; IQ<=18; IQ++){
              if(eN[MODE][32][IQ]>0) {
                  fprintf(f09,"%5d|",eN[MODE][32][IQ]);
                  s << "<td>" << QString::number(eN[MODE][32][IQ]) << "</td>";
                }else    {                fprintf(f09,"     |");
                  s << "<td></td>";
                }
            }
          s << "</tr>";
          fprintf(f09,"\\par | dSig/dOmega  | ");
          s << "<tr><td>dSig/dOmega</td>";
          for(I=1; I<=18; I++) {
              TET=(I-0.5)*0.1745;
              _DSIG[I]=eN[MODE][32][I]/(_NCASC*sin(TET)*0.1745*6.2832)*_SIGMA;
              if(_DSIG[I]>0) {
                  fprintf(f09,"%5.0g|",_DSIG[I]);
                  s << "<td>" << QString::number(_DSIG[I]) << "</td>";
                }else    {
                  fprintf(f09,"  0  |");
                  s << "<td>  0 </td>";
                }
            }
          s << "</tr>";
          fprintf(f09,"%s",s500);
          // s << QString::fromUtf8(s500html);
          if(IWE==0)continue;


          for(KW=1; KW<=3; KW++)    {
              fprintf(f09,"\\par |%5.1f - %5.1f | ",eEW1[KW],eEW2[KW]);
              s << "<tr><td>" << QString::number(eEW1[KW]) << " - " << QString::number(eEW2[KW]) << "</td>";
              for(IQ=1; IQ<=18; IQ++)
                if(eNW[MODE][KW][IQ]>0) {
                    fprintf(f09,"%5d|",eNW[MODE][KW][IQ]);
                    s << "<td>" << QString::number(eNW[MODE][KW][IQ]) << "</td>";
                  }else {
                    s << "<td></td>";
                    fprintf(f09,"     |");
                  }
              s << "</tr>";
            }

          s << " <tr><td>Above " << QString::number(eEW1[4]) << " </td> ";
          fprintf(f09,"\\par | Above %6.1f | ",eEW1[4]);

          for(IQ=1; IQ<=18; IQ++){
              if(eNW[MODE][4][IQ])  {
                  fprintf(f09,"%5d|",eNW[MODE][4][IQ]);
                  s << "<td>" << QString::number(eNW[MODE][4][IQ]) << "</td>";
                } else  {
                  s << "<td></td>";
                  fprintf(f09,"     |");
                }
            }
          s << "</tr>";
          //s << "</p><p>" << QString::fromUtf8(s500html) << "</p>";
          fprintf(f09,"%s\\fs16 ",s500);
          s << "</table>";
        }


      //---------------------------------------- L36


      MODE=4;
      fprintf(f09,"\\par\\par\\b\\f2\\cf1\\fs16 %s     Gamma\\cf3  ray spectrum \\cf0 %d\\cf1  events\\b0\\f2\\cf2\\fs16 ",s35x,eNT[4]);
      s << "<p>&nbsp;</p><p><span style=\"color:red\">  Gamma </span>ray spectrum " << QString::number(eNT[4]) << "  events</p>";


      fprintf(f09,"\\par\\par  Emission from unbound and bound states(*), and total gamma ray spectrum");
      fprintf(f09,"\\par  (*) note that emission of a particle from an unbound state is not allowed"
                  " in the code if Ecm is less than Emin \\f0\\fs14\\cf0\\par ");
      s << "<p>&nbsp;</p><p> Emission from unbound and bound states(*), and total gamma ray spectrum</p> ";
      s << "<p>  (*) note that emission of a particle from an unbound state is not allowed"
           " in the code if Ecm is less than Emin </p>";
      // fprintf(f09,"%s",s500a);
      s << QString::fromUtf8(s500ahtml);
      fprintf(f09,"\\par |        \\b    Energy range (MeV)  \\b0      |  \\b    Unbound  \\b0     |   \\b     Bound   \\b0 |     TOTAL   |"
                  "%s", s502);
      s << "<table border=1 cellspacing=0 cellpadding=2><tr><th>  Energy range (MeV) </th><th> Unbound </th><th> Bound </th><th> TOTAL </th></tr>";
      for(K=1; K<=30; K++) {
          E1=K-1;
          E2=E1+1.;
          if(eNGT[K]>0){
              fprintf(f09,"\\par |             %6.1f - %-6.1f        |     %6d       |     %6d    |    %6d   |",
                      E1,E2,eNG[K],eNGB[K],eNGT[K]);
              s << "<tr><td align=center>   " << QString::number(E1) << " - " << QString::number(E2) << " </td><td> " << QString::number(eNG[K]) << " </td><td> " << QString::number(eNGB[K]) << "</td><td>" << QString::number(eNGT[K]) << " </td></tr>";
            }
        }


      E1=30.*DELE;
      if(eNGT[31]>0){
          fprintf(f09,"\\par |            Above  %6.1f            |     %6d       |     %6d    |    %6d   |",
                  E1,eNG[31],eNGB[31],eNGT[31]);
          s << "<tr><td align=center> Above " << QString::number(E1) << " </td><td> " << QString::number(eNG[31]) <<  "</td><td> " << QString::number(eNGB[31]) << "</td><td> " << QString::number(eNGT[31]) << "</td><td></tr>";
        }
      fprintf(f09, "%s", s502);
      // s << QString::fromUtf8(s502html);
      //E1=0.;

      fprintf(f09,  "\\par |    \\b         Total\\b0                    |     %6d       |     %6d    |    %6d   |",
              eNG[32],eNGB[32],eNGT[32]);
      s <<"<tr> <th>  Total </th><td>  " << QString::number(eNG[32]) << "</td><td> " << QString::number(eNGB[32]) <<  "</td><td>  " << QString::number(eNGT[32]) << "</td></tr>";
      fprintf(f09,"%s",s500a);
      //s << QString::fromUtf8(s500ahtml);
      s << "</table>";
      // s.flush();
    }
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWw
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::OPTPOT(int I,int IN, int IZ)
{
  //     WRITTEN BY A. GAVRON
  //=========      common/POT/
  extern double _V0D[5],_R0RD[5],_ARD[5],_R0CD[5],/*_RCCC[4],_VQS[4],_VQ0[4],*/_W0D[5],_R0ID[5],_AID[5],
      _RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
  extern int _NPD[5],_IMAG[5],_IRAD[5];

  //   - SUBROUTINE TO PROVIDE DEFAULT OPTICAL MODEL PARAMETERS
  double A=IN+IZ;

  switch (I)  {
    case 1 :
      //   - NEUTRONS
      _V0D[1]=47.01;
      _V01D[1]=-.267;
      _V02D[1]=-.0018;
      //
      _R0RD[1]=1.322-A*7.6E-4+pow(A,2)*4.E-6-pow(A,3)*8.E-9;
      _ARD[1]=.66;
      _R0CD[1]=0.;
      _W0D[1]=9.52;
      _W01D[1]=-.053;
      _W02D[1]=0.;
      _R0ID[1]=1.266-A*3.7E-4+pow(A,2)*4.E-6-pow(A,3)*4.E-9;
      _AID[1]=.48;
      _RMCHD[1]=0.;
      _NPD[1]=250;
      _IMAG[1]=1;
      _IRAD[1]=1;
      break;

    case 2 :

      //     PROTONS
      //     PROTONS AND NEUTRONS FROM PEREY AND PEREY
      _V0D[2]=53.3+27*(IN-IZ)/A+.4*IZ/pow(A,.3333);
      _V01D[2]=-.55;
      _V02D[2]=0.;
      _R0RD[2]=1.25;
      _ARD[2]=.65;
      _R0CD[2]=1.25;
      _W0D[2]=13.5;
      _W01D[2]=0.;
      _W02D[2]=0.;
      _R0ID[2]=1.25;
      _AID[2]=0.47;
      _RMCHD[2]=0.;
      _NPD[2]=250;
      _IMAG[2]=1;
      _IRAD[2]=1;
      break;

    case 3 :
      //     ALPHAS FROM IGO AND HUIZENGA
      _V0D[3]=50.;
      _V01D[3]=0.;
      _V02D[3]=0.;
      _R0RD[3]=1.17*pow(A,.3333)+1.77;
      _ARD[3]=0.576;
      _R0CD[3]=_R0RD[3]-1.77;
      _W0D[3]=3.0+.105*A;
      _W01D[3]=0.;
      _W02D[3]=0.;
      _R0ID[3]=_R0RD[3];
      _AID[3]=_ARD[3];
      _RMCHD[3]=0.;
      _NPD[3]=250;
      _IMAG[3]=0;
      _IRAD[3]=0;
      break;
    };
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void TLLL(double ER, int NP,int IZ,int IN, double TL[25], int &LLMAX, QTextStream& s)
{
  //     WRITTEN BY A. GAVRON
  //=======================================================
  //=========      common/MASS/
  //  extern double _FLMASS,_ARATIO;
  //   extern double _FBARR[111];
  //extern int _NFILE;  ///////////////MICHELLE
  //   extern int _LFISS;
  //   extern int _IZPART[5], _INPART[5];
  //=========      common/POT/
  extern double /*_V0D[5],_R0RD[5],_ARD[5],_R0CD[5],*/_RCCC[4],_VQS[4],_VQ0[4];//_W0D[5],_R0ID[5],_AID[5],
  //   _RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
  //   extern int _NPD[5],_IMAG[5],_IRAD[5];
  //=========      common/SRCH/
  //  extern double _PROB[5556],_EMAX[5],_EMIN[5],_BE[5],_GE1,_GE2;
  extern double _EMAX[5],_EMIN[5];//,_BE[5],_GE1,_GE2;
  //    extern int    _IPROB,_MJ[5556],_IMODL;
  //=========      common/OLEG/
  extern int /*_INPUT,*/ _IDIST;
  //==========================================================
  // _IDIS is responsible for PRINT - analog CTLL

  extern void TLOM(int IPOT, double AMP, double AZP, double AMT,
                   double AZT, double ECM, int &LL, double TL[Max_MOM+1]);

  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  static double  TLPN[3][13][31],TLAL[25][31],TLT[Max_MOM+1],DE[4],LL[4][31], DE1[4];

  const char *NMODE[5][7]={"    ","Neutron","Proton","Alpha","HI"};
  char s52[128]="\\par\\par\\b  Mode = \\cf3 %s\\cf0\\b0 "
                "\\par\\b   EP / L       0     1     2     3     4     5     6     7     8     9    10    11\\b0 ";
  char s52html1[53] = "<p>&nbsp;</p><p><b> Mode = <span style=\"color:blue\">";
  char s52html2[169] = "</span></b></p><table><tr><th>EP / L</th><th>0</th><th>1</th><th>2</th><th>3</th><th>4</th><th>5</th><th>6</th><th>7</th><th>8</th><th>9</th><th>10</th><th>11</th></tr>";
  //     CHANGE 'CTLL' TO BLANKS TO OBTAIN TRANSMISSION COEFF PRINTOUT

  //      double PROB;
  double TLL=-87.49;      // check K for increasing
  double ZP=0, AP=0, E=0, EE=0, DEE=0, DTT=0, T0=0;
  // int LC, IP, IEN, L, NE, LL1, LP1, LMT, KP, LC1=0;
  int LC=0, IP=0, L=0, NE=0, LL1=0, LP1=0, LMT=0, KP=0, LC1=0;

  //=========================================================================
  if(_cl_K<=0) {
      if(_IDIST==2) {
          fprintf(f09,"\\par\\par\\b      Transmission coefficients for particle emission\\b0 ");
          s << "<p>&nbsp;</p><p><b>  Transmission coefficients for particle emission</b></p> ";
        }
      double A=IZ+IN;
      double Z=IZ;
      _cl_K=1;
      bool tab = false;
      for(IP=1; IP<=2; IP++) {

          DE[IP]=(_EMAX[IP]-_EMIN[IP])/29.;
          ZP=0.;
          AP=1.;
          if(IP==2)ZP=1.;
          if(_IDIST==2){
              tab = true;
              fprintf(f09,s52,NMODE[IP]);
              s << s52html1 << NMODE[IP]<< s52html2;
            }
          for(int IEN=1; IEN<=30; IEN++) {

              E=_EMIN[IP]+(IEN-1)*DE[IP];
              TLOM(IP,AP,ZP,A-1.,Z-1.,E,LC,TLT);
              LC1=LC+1;
              if(_IDIST==2) {
                  fprintf(f09,"\\par  %6.2f    ",E);
                  s << "<tr><td>  " << QString::number(E) << "</td>";
                  for(KP=1; KP<=LC1; KP++) {
                      if(KP!=1 && ( (KP-1)%12 ==0 ) ){
                          s << "</tr><tr><td></td>   ";
                          fprintf(f09,"\\par            ");
                        }
                      fprintf(f09,"%6.3f",TLT[KP]);
                      s << "<td>" << QString::number(TLT[KP]) << "</td>";
                    }
                  s << "</tr>";
                }
              LL[IP][IEN]=11;
              if(LC<11)LL[IP][IEN]=LC;
              for(int L=1; L<=12; L++) {
                  TLPN[IP][L][IEN]=-87.49;
                  if(L<=LC+1 && TLT[L]>0.)TLPN[IP][L][IEN]=log(TLT[L]+1.E-38);
                }
            }
          //if(_IDIST == 2) {s << "</table>";}
          if(tab == true) { s << "</table>"; }
        }
      tab = false;
      AP=4.;
      ZP=2.;
      DE[3]=(_EMAX[3]-_EMIN[3])/29.;
      IP=3;
      if(_IDIST==2){
          tab = true;
          fprintf(f09,s52,NMODE[IP]);
          s << s52html1 << NMODE[IP] << s52html2;
        }

      for(int IEN=1; IEN<=30; IEN++) {
          E=_EMIN[3]+(IEN-1)*DE[3];
          TLOM(3,AP,ZP,A-4.,Z-2.,E,LC,TLT);
          LC1=LC+1;
          if(_IDIST==2) {
              fprintf(f09,"\\par  %6.2f    ",E);
              s << "<tr><td> " << QString::number(E) << "</td>";
              for(KP=1; KP<=LC1; KP++) {
                  if(KP!=1 && ( (KP-1)%12 ==0 ) ) {
                      fprintf(f09,"\\par            ");
                      s << "<tr><td></td>    ";
                    }
                  fprintf(f09,"%6.3f",TLT[KP]);
                  s << "<td>" << QString::number(TLT[KP]) << "</td>";
                }
              s << "</tr>";
            }
          LL[3][IEN]=LC;
          if(LC>23)LL[3][IEN]=23;
          for(L=1; L<=24; L++) {
              TLAL[L][IEN]=log(1.E-38);
              if(L<=LC+1 && TLT[L]>0.) TLAL[L][IEN]=log(TLT[L]+1.E-38);
            }
        }
      if(tab == true){s << "</table>";}
      for(IP=1; IP<=3; IP++) DE1[IP]=1./DE[IP];

      _VQ0[1]=0.;
      _VQ0[2]=1.44*(Z-1.)/(_RCCC[2]*pow((A-1.),0.3333));
      _VQ0[3]=2.88*(Z-2.)/(_RCCC[3]*pow((A-4.),0.3333));

    }
  //=========================================================================
  else  {
      EE=ER+_VQ0[NP]-_VQS[NP];
      //     THIS PROCEEDURE ADJUSTS THE T-LS TO THE CHANGING COULOMB
      //     BARRIER DURING DEEXCITATION.
      if(EE>=_EMAX[NP])EE=_EMAX[NP]-.01;
      if(EE<=_EMIN[NP])EE=_EMIN[NP]+.01;

      NE=(EE-_EMIN[NP])*DE1[NP];
      DEE=EE-_EMIN[NP]-DE[NP]*NE;
      NE++;
      LLMAX=LL[NP][NE];
      LL1=LLMAX+1;
      if(NP==3)  {
          TL[1]=0.;
          LMT=LLMAX+1;
          for(LP1=1; LP1<=LMT; LP1++) {
              TL[LP1]=0.;
              T0=TLAL[LP1][NE];
              DTT=(TLAL[LP1][NE+1]-T0)*DE1[NP];
              TLL=T0+DTT* DEE;
              TL[LP1]=exp(TLL);
            }
          _cl_K++;
        }
      else
        for(LP1=1; LP1<=LL1; LP1++) {
            TL[LP1]=0.;
            T0=TLPN[NP][LP1][NE];
            DTT=(TLPN[NP][LP1][NE+1]-T0)*DE1[NP];
            TLL=T0+DTT* DEE;
            TL[LP1]=exp(TLL);
          }

    }
  //=========================================================================
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void PACE::MOMENT(int ID,double A, double AP, int MODE, double EP,
                  fusion_event *f, QTextStream& s)
{                                           //      EP -- particles energy in CM
  //
  //     ORIGINAL JULIAN SUBROUTINE
  //
  //     ***********
  //      MONITORING OF OUTGOING PARTICLES.
  //      common/MOM/
  //=========      common/XQANG/
  extern double _CT;       // COSTH  (fortran source)
  extern int    _MDIR;
  //-------------------------------------

  if(ID < 2 ) {                           //   INITIALIZATION
      f->y=EP*_MDIR;
      f->x=0.;
      f->z=EP*(1.-_MDIR);
    }
  else if(ID==2)  {                       //   CALCULATION OF MOMENT

      double RN4=RANF();
      double PHI=6.28318*RN4;
      double SOX=2.*AP*EP/(A*(AP+A)*938.);
      //oleg      if(SOX<0.) {
      //oleg          fprintf(f09,"\\parAP,EP,A =  %12.4e %12.4e %12.4e",AP,EP,A);
      //    oleg      return; // oleg
      //oleg          }

      double VT=mzsqrt(SOX);
      SOX = max(1. - _CT*_CT, 0);
      double SINTH=mzsqrt(SOX);

      double VZSE=VT*_CT;
      double VYSE=VT*SINTH*sin(PHI);
      double VXSE=VT*SINTH*cos(PHI);

      f->z+=VZSE;
      f->y+=VYSE;
      f->x+=VXSE;

      double VZP=f->z - VZSE*(A/AP+1.);
      double VYP=f->y - VYSE*(A/AP+1.);
      double VXP=f->x - VXSE*(A/AP+1.);

      double VPP2 = VZP*VZP + VYP*VYP + VXP*VXP;

      f->particle_energy=AP*VPP2*469.;      // LAB
      SOX=0.;
      if (VPP2 > 0.0) SOX = VZP/ sqrt(VPP2);  //     DO NOT USE QUICK FUNCTIONS HERE
      if (fabs(SOX)>1.) {
          fprintf(f09,"\\par  SOX= %f",SOX);
          s << "<p>  SOX = " << QString::number(SOX) << "</p>";
        }
      f->particle_angle = acos(SOX)*GRAD;
      //s.flush();
      OUTEM(2,s,MODE,f->particle_energy,f->particle_angle);
    }
  else if(ID>2) {   //======================================================================


      double VFTS=f->z*f->z+f->x*f->x+f->y*f->y;             // oleg - check this equivalence
      f->y=0.5*A*VFTS*938.;
      if(VFTS==0.)     f->z=0;
      else             f->z= acos(f->z/ sqrt(VFTS))*GRAD;
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::TRACK(int K,int JI,int JF,double EI,double EF,double ENERGY,
                 int IC,double SPROB,double ECLOSS, QTextStream& s)
{                                                           // ic = 0 - make 0
  //=========      common/OUT/                                // ic = 1 - init E
  extern double _RLEV[4][DIM_RLEV];       // TRANSPORTED     // ic = 2 - monitoring
  // ic = 3 - print
  //=========      common/SCH/
  extern double _EBIN[5][DIM_EBIN];         // TRANSPORTED
  extern int    _MAXJ[5][DIM_EBIN],_MAXJS[5][DIM_EBIN];   // TRANSPORTED
  //=========      common/QEJS/
  extern int _MEBIN[5];


  //     WRITTEN BY A. GAVRON
  static int  NP[6][4],JD[6][4],JD2[6][4],NTOT[4],JAVR[4],JRMS[4];
  static double DE[6][4],GAMMV[4],TAUV[4];

  double  ATOT, AJAVR, AJRMS, AN, A, B, C;
  static double E[4];

  const char *PART[6][6]={"      ","Neut","Prot","Alph","Gamm","Fiss"};
  QString part_qt[6] = {"      ","Neut","Prot","Alph","Gamm","Fiss"};
  //double RLEV;
  int I,KI,J=0, MEB=0, MI=0,INDX=0;
  double GAMMA,TAU;

  //======================================================
  //======================================================
  switch(IC) {

    //======================================================
    case 0 :
      for(I=0; I<4; I++) {
          E[I]=NTOT[I]=JAVR[I]=JRMS[I]=GAMMV[I]=TAUV[I]=0;
          for(KI=0; KI<6; KI++) {
              // NP[K][KI]=JD[K][KI]=JD2[K][KI]=DE[K][KI]=0;
              NP[KI][K] =0;JD[KI][K]=0;JD2[KI][K]=DE[KI][K]=0;
            }
        }
      break;
      //======================================================
    case 1 :
      E[1]=ENERGY-ECLOSS-0.1;
      E[2]=ENERGY*0.5+0.5;
      E[3]=10.0;
      break;

      //======================================================
    case 2 :
      J=0;
      if(EI>E[1]&&EF<=E[1])      J=1;
      else if(EI>E[2]&&EF<=E[2])      J=2;
      else if(EI>E[3]&&EF<=E[3])      J=3;

      if(J*SPROB==0)return;

      MEB=_MEBIN[4];

      for(I=1; I<=MEB; I++) {
          MI=MEB-I+1;
          if(fabs(_EBIN[4][MI]-ENERGY) <0.501) break;       // TRANSPORTED
        }

      if(I>MEB){
          fprintf(f09,"\\par Track error. LD range = %10.4f  %10.4f    Energy = %10.4f   J = %d",
                  _EBIN[4][1],_EBIN[4][MEB],ENERGY,J);     // TRANSPORTED
          s << "<p> Track error. LD range = " << QString::number(_EBIN[4][1]) << " " << QString::number(_EBIN[4][MEB]) << "  Energy = "
                                                                                                                       << QString::number(ENERGY) << "  J = " << QString::number(J);
        }

      INDX=_MAXJS[4][MI]-_MAXJ[4][MI]+JI;                   // TRANSPORTED
      GAMMA=SPROB/(_RLEV[3][INDX]*6.2832+1.E-33)+1.E-33;   // TRANSPORTED     RLEV[4] it was before
      TAU=6.6E-22/GAMMA;

      GAMMV[J]+=GAMMA;
      TAUV[J]+=TAU;
      NTOT[J]++;
      JAVR[J]+=JI;
      JRMS[J]+=JI*JI;
      NP[K][J]++;
      JD[K][J]+=JF-JI;
      JD2[K][J]+=(JF-JI)*(JF-JI);
      DE[K][J]+=EI-EF;

      break;
      //===============================================================================

    case 3 :
      fprintf(f09,"\\par%s\\par  Track down of decay modes at \\b %.0f, %.0f,\\b0  and\\b  %.0f\\b0  MeV excitation",
              s40star, E[1],E[2],E[3]);
      s << "<p> Track down of decay modes at <b> " << QString::number(E[1]) << ", " << QString::number(E[2]) << ", " << QString::number(E[3]) << "</b>  MeV excitation </p>";
      for(I=1; I<=3; I++) {
          ATOT=NTOT[I]+1.E-9;
          GAMMV[I]/=ATOT;
          TAUV[I]/=ATOT;
          AJAVR=JAVR[I];
          AJRMS=JRMS[I];
          AJAVR/=ATOT;
          AJRMS=sqrt(AJRMS/ATOT-AJAVR*AJAVR);
          fprintf(f09,"\\par\\par\\b  Ex = %.0f\\b0\\fs14   Gamma = %.2e MeV   Lifetime = %.2e sec"
                      "   Average J = %.1f   Stand.dev. = %.1f\\fs16 "
                      "\\par\\par\\cf2               Part    Num   DelJ    RMS-dJ\\cf0 ",
                  E[I],GAMMV[I],TAUV[I],AJAVR,AJRMS);
          s << "<p>&nbsp;</p><table cellpadding='15'><tr><th> Ex = " << QString::number(E[I]) << " </th><th> Gamma = " << QString::number(GAMMV[I],'g',3) << " MeV </th><th>   Lifetime = "
            << QString::number(TAUV[I],'g',3) << " sec  </th><th>Average J = " << QString::number(AJAVR,'f',3) << " </th><th>Stand.dev. = " << QString::number(AJRMS, 'f',3)
            << "</th></tr></table><p>&nbsp;</p><table cellpadding='8'><tr><th></th><th>Part</th><th>Num</th><th>DelJ</th><th>RMS-dJ</th></tr>";
          for(K=1; K<=5; K++) {
              AN=NP[K][I];
              if(AN<=0.)continue;
              A=JD [K][I]/(AN+1.E-9);
              // B=fabs(JD2[K][I]);
              B=abs(JD2[K][I]);
              B=sqrt(B/(AN+1.E-9));
              C=DE[K][I]/(AN+1.E-9);
              if(K <5){
                  fprintf(f09,"\\par       %s  %6.0f%7.1f %7.1f %8.1f",(const char*)PART[K],AN,A,B,C);
                  s << "<tr><td>      " << part_qt[K] << "</td><td align=center> " << QString::number(AN) << "</td><td align=center>" << QString::number(A)
                    << "</td><td align=center>  " << QString::number(B) << "</td><td align=center>" << QString::number(C) << "</td></tr>";
                }else{
                  fprintf(f09,"\\par       %s  %6.0f",(const char*)PART[K],AN);
                  s << "<tr><td>     " << part_qt[K] << "</td><td>   " << QString::number(AN) << "</td></tr>";
                }

            }
          s << "</table>";
        }
      //===============================================================================
    }

  //===============================================================================
}
//===============================================================================


