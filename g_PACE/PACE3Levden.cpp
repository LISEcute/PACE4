#include <QTextStream>
#include <QString>
//#include "WinLise_Constant.h"
#include "ftype.h"
#include <math.h>
#include <stdio.h>
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#include "ParticleStates.h"
#include "pace.h"

extern  int _cl_IPRNT;
extern  int _cl_KK;
extern QString qbuf;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
extern double RANF();
extern void TLLL(double ER, int NP,int IZ,int IN, double TL[25], int &LLMAX, QTextStream& s);
extern double MASSES(int IZ, int IN, int &opt1);
//extern double SIERK(double A, double Z, double AL, double BARFAC);

extern FILE *f09;
//extern QTextStream htmlStream;

extern const char *s40star;
extern const char *s10star;
//extern char *s35x;
const char s35x[10] = "         ";
extern char *s100m;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//=========      common/MASS/
extern double _FLMASS,_BARFAC,_ARATIO,_BARMAX,_FBARR[Max_MOM+1];
extern int    _IZR[5],_INR[5];
//_NFILE,
extern int _LFISS, _IZPART[5], _INPART[5];

//=========      common/QMMJ/
extern double _FYRST,_FACLA,_EROT[Max_MOM+1];
extern int    _MAXC;

//=========      common/DAT/
extern double _SPART[5], _PMASS[5];
//=========      common/GAM/
extern double  _FGE1r,_FGM1r,_FGE2r,_FGM2r;
extern double _FGE1,_FGM1,_FGE2,_FGM2,_DISCPR,_ERGMIN ,_EDL[37][7];
extern int  _IADEF,_NDISC,_NDL[7],_IZDL[7],_INDL[7],_JDL[37][7];
//=========      common/XQRES/
extern double _ENERGY,_VZC;
extern int    _IZ,_IN,_ITRAC;
//=========      common/SRCH/
extern double _PROB[5556],_EMAX[5],_EMIN[5],_BE[5],_GE1,_GE2,_GM1,_GM2;
extern int    _IPROB,_MJ[5556],_IMODL;
//==================     common /XQTT/
extern double _UX,_EX,_TT;
//-------------  common/RLVV/
double RLEVT[DIM_RLEV];

//=========      common/SCH/
extern double _EBIN[5][DIM_EBIN],_RMASS[5];      // OLEG
extern int    _MAXJ[5][DIM_EBIN],_MAXJS[5][DIM_EBIN]; // OLEG
//=========      common/QEJS/
extern int _MEBIN[5];
//=========      common/XQANG/
extern double _SUM[11][301],_CT;
extern int    _MDIR;
//=========      common/FAC/
extern double FA[204];

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void AMASS(double *A, fusion_event &c, int MODE, int &MEBIN, int *MAXJ, int *MAXJS,
                 double *EBIN, double *RLEV, int &IXPR, double *BE, QTextStream &s);
void RANGE(int MODE,int &IMIN,int &IMAX, double Ex);

void FIND(double E, int &IXRSM, double *EBIN, int MEBIN, int MODE, QTextStream &s);


void GCLVD(double &ALIT, double /*SIGSQ*/, double FACT1, double DELTA, double EFIRST,
     double EMAX, double DELEX, double SR, int &KK, int &MEBIN, int IXPR, int *MAXJ,
     int *MAXJS, double *EBIN, double *RLEV, QTextStream &s);
void SEARCH(fusion_event *a, fusion_event &c, int &MODE, QTextStream &s);
void MJRAN(fusion_event &csi, fusion_event &csf, int MODE, double EP, QTextStream &s);
double C3J(double FJ1,double FJ2,double FJ3,double FM1,double FM2,double FM3);


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::AMASS(double *A, fusion_event &c, int MODE, int& MEBIN, int *MAXJ, int *MAXJS,
                 double *EBIN, double *RLEV, int &IXPR, double *BE, QTextStream& s)
  {
//
//     MODIFIED JULIAN SUBROUTINE
//
//     **************
//      DETERMINE MASSES.

//======================================================================================

      int MAXXL=_MAXC+9;
      MAXXL = min(110,MAXXL);

      _IZR[MODE]=c.Z-_IZPART[MODE];
      _INR[MODE]=c.N-_INPART[MODE];

      int opt1;
      A[MODE] = MASSES(_IZR[MODE],_INR[MODE], opt1);

      BE[MODE]=A[4]-A[MODE]-_PMASS[MODE];

      if(MODE==4) {
              double AA=_IZR[MODE]+_INR[MODE];
              double ZZ=_IZR[MODE];

              for(int I=1; I<=MAXXL; I++)
                      _FBARR[I]=SIERK(AA,ZZ,I-1,_BARFAC);
              }


      int IMIN=0; int IMAX=0;

      RANGE(MODE,IMIN,IMAX, c.Ex);
      LEVDEN(_IZR[MODE],_INR[MODE],IMIN,IMAX,MAXJ,MAXJS,EBIN,MEBIN,RLEV,IXPR,BE[MODE],1.,s);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void RANGE(int MODE,int &IMIN,int &IMAX, double Ex) {

//     WRITTEN BY A. GAVRON
//   COMMON/SRCH/PROB(5555),IPROB,MJ(5555),QMAX(4),QMIN(4),BE(4),IMODL
//=========      common/SRCH/
//extern double _PROB[5556],_EMAX[5],_EMIN[5],_BE[5],_GE1,_GE2;
//extern int    _IPROB,_MJ[5556],_IMODL;

      double EMAX;
      double EMIN, EMIN1, BR;


      // check OLEG       int IZR=_IZR[MODE];
      // check OLEG       int INR=_INR[MODE];
      // check OLEG       double A=IZR+INR;
      // check OLEG       double Z=IZR;

      int IMM=Max_EBIN;
      int INN=IMM-1;

   EMIN = _ERGMIN + _BE[MODE] - _EMAX[MODE] - 2.*_EMAX[4];
   EMAX = Ex + _BE[MODE] - _EMIN[MODE];

if(MODE==4) {
                //     THE RANGE FOR MODE=4 ALLOWS FOR 4 GAMMA DECAYS. IF YOU SPECIFY
                //     A LARGE GAMMA DECAY WIDTH, FIND-ERRORS MAY RESULT. INCREASE THE
                //     COEFFICIENT OF QMAX(4) TO INCREASE LEVEL DENSITY RANGE.

      EMIN-=_EMAX[4];
      BR=max(_FBARR[1],_BARMAX)+1.;
      EMIN1=_ERGMIN-BR-15.;

                //     THIS ASSURES ENOUGH LEVEL DENSITY FOR FISSION CALCULATION
      EMIN = min(EMIN1,EMIN);
      }

      IMIN = EMIN-2.5;
      IMIN = max(0,IMIN);
      IMAX = EMAX+2.5;

      IMAX = min(IMAX,IMM);
      IMIN = min(IMIN,INN);

      if(IMAX<=0) IMAX=2;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::LEVDEN(int IZ, int IN, int IFIRST, int IMAX, int *MAXJ, int *MAXJS, double *EBIN,
 int &MEBIN, double *RLEV, int &IXPR, double /*BE*/, double ARATIO, QTextStream &s)
{
//
//     ORIGINAL JULIAN SUBROUTINE
//
//     ***********
//      GILBERT AND CAMERON DEFAULT LEVEL DENSITY SUBROUTINE

//      COMMON/GAM/FGE1,FGM1,FGE2,FGM2,ERGC,DISCPR,IADEF,ND,ERGMIN        PAC41200
//     X ,NDL(6),IZDL(6),INDL(6),JDL(36,6),EDL(36,6)                      PAC41210
//      COMMON/QMMJ/MAXC,FYRST,FACLA,EROT(110)                            PAC41220


int I, J, MAX;
int JBIN[DIM_EBIN] = {0};
int IA=IZ+IN;

if(IA*IZ<=0)
     return;

int KK=0;
int IBIN=1;
IFIRST=max(IFIRST,1);
IXPR=0;


if(_cl_IPRNT==0) {
      fprintf(f09,"\\par\\par\\cf2\\f2\\fs20  ");
      s << "<p>&nbsp;</p>";
      if(_FYRST==0.){
          fprintf(f09," *** Gilbert - Cameron spin cutoff parameter used");
          s << "<p> *** Gilbert - Cameron spin cutoff parameter used</p>";
      }else {
          fprintf(f09," *** C.P.S. rotating liquid drop yrast line used \n multiplied by factor of %.2f",_FYRST);
          s <<  "<p>*** C.P.S. rotating liquid drop yrast line used <br> multiplied by factor of " << qbuf.number(_FYRST) << "</p>";
      }
          fprintf(f09,"\\cf0\\f0\\fs16  ");

      _cl_IPRNT=1;
      }


if(IFIRST<=1) {
      IXPR=1;
      JBIN[1]=1;
      EBIN[1]=0;
      if(_NDISC!=0)
              for(I=1; I<=_NDISC; I++) {
                      if(IZ!=_IZDL[I] || IN!=_INDL[I])continue;
                      IXPR=_NDL[I];
                      for(J=1; J<=IXPR; J++) {
                              JBIN[J]=_JDL[J][I];
                              EBIN[J]=_EDL[J][I];
                              }
                      break;
                      }

      }


//L1:   format(10I5);
//L4:   format((8(F6.0,1X,I3)));

if(IXPR!=0)
   for(I=1; I<=IXPR; I++) {
      MAX=JBIN[I];
      MAXJ[I]=MAX;
      for(J=1; J<=MAX; J++) {
              KK++;
              RLEV[KK]=0.;
              if(J==MAX)RLEV[KK]=1.E-38;
              }
      MAXJS[I]=KK;
      }

      double AMR=IA;
      double AZR=IZ;
      double DELEX=IBIN;
      double DEL2=DELEX*0.5;
      if(IXPR>0)IFIRST=(EBIN[IXPR]+0.5+DEL2*0.5);
      double EFIRST=DEL2+IFIRST;
      double EMAX=IMAX;
      double SR=((IA)%(2))*0.5;

      double FALIT=0, ALIT=0, SIGSQ=0, FACT1=0, DELTA=0, ACONST=0;

    if(AMR*AZR > 0 ){
       GCLDP(AMR,AZR,_IADEF,FALIT,ALIT,SIGSQ,FACT1,DELTA,ACONST,SR, ARATIO, s);
       GCLVD(ALIT,SIGSQ,FACT1,DELTA,EFIRST,EMAX,DELEX,SR,KK,MEBIN,IXPR,MAXJ,MAXJS,EBIN,RLEV,s);
       if(ALIT<=0.)  {
           fprintf(f09,"\\par ALIT<=0 Z=%d A=%d IFIRST=%d IMAX=%d\\par",IZ,IA,IFIRST,IMAX);
           s << "<p> ALIT<=0 Z = " << qbuf.number(IZ) << " A = " << qbuf.number(IA) << " IFIRST = " << qbuf.number(IFIRST)
                << " IMAX = " << qbuf.number(IMAX) << "</p>";
       }
    }
    else {
        fprintf(f09,"\\par LEVDEN: A*Z <= 0 \\par");
        s << "<p> LEVDEN: A*Z <= 0  </p>";
    }



//OLEG      else          MAXXL=_MAXC+9;



//     WRITE(9, 82)IZ,IA,IXPR,IFIRST,IMAX,ALIT,MAXXL,EROT(MAXXL),DELTA,BE
//L82:  format(5I5,F7.1,I5,F7.2,F6.2,F9.3);
//     IF(IXPR.EQ.0)RETURN
//     WRITE(9, 494)(EBIN(I),JBIN(I),I=1,IXPR)
//     WRITE(9, EXPERIMENTAL LEVELS (IF EXIST)
// 494 FORMAT((8(1X,F8.3,1X,I3)))
s.flush();
      }

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void GCLVD(double &ALIT, double /*SIGSQ*/, double FACT1, double DELTA, double EFIRST,
     double EMAX, double DELEX, double SR, int &KK, int &MEBIN, int IXPR, int *MAXJ,
     int *MAXJS, double *EBIN, double *RLEV, QTextStream& s)
 {
//     ***********
      double FACROT[111],TWOSR[111];


//      common/QMMJ/MAXC,FYRST,FACLA,EROT[110];
//      common /XQTT/UX,EX,TT;
//      common/GAM/FGE1,FGM1,FGE2,FGM2,ERGC,DISCPR,IADEF;

      int J,IEX, IE, JMX, JJ;
      int II=Max_EBIN-IXPR;
      int MAXXL=min(_MAXC+9,110);
      double AJR,EEXRS,EEXRST, UEX, TEMP, FACT2, FACTT, AVLOG,RLEVD1;


      MEBIN=Max_EBIN;

for(J=1; J<=MAXXL; J++) {
              AJR=J-1.+SR;
              FACROT[J]=1.;
              if( !(_IADEF==2 || AJR>22) ) {
                        if(AJR>=18.) FACROT[J]=0.5+(AJR-18.)*0.125;
                        else         FACROT[J]=0.5;
                        }

              FACROT[J]=1./FACROT[J] *_EROT[J];
              }

//--------------------------cikl start L6
for(IEX=1; IEX<=II; IEX++) {
      IE=IEX+IXPR;
      EEXRS=EFIRST+(IEX-1)*DELEX;
      EEXRST=EEXRS;
      if(EEXRST<_EX)EEXRS=_EX;

      if(EEXRST>EMAX) {MEBIN=IE-1; break;}

      EBIN[IE]=EEXRST;
      UEX=EEXRS-DELTA;
      if(UEX<=0 ) {
              MAXJ[IE]=1;
              KK++;
              RLEV[KK]=0.;
              MAXJS[IE]=KK;
              continue;
              }

      if(ALIT<=0)
                {ALIT=0; return;}          // OLEG 06/13/2011

      TEMP=sqrt(UEX/ALIT);

      FACT2=FACT1/(pow(UEX,1.25)*pow(TEMP,1.5));
      if(EEXRST<=_EX) FACTT=exp((EEXRST-_EX)/_TT);
      else            FACTT=1.;


      FACT2*=FACTT;
      MAXXL = min(MAXXL,110);
                        //    FOLLOWING LOOP REPLACED FOR OPTIMIZATION ENHANCEMENT
      JMX=MAXXL;

      for(J=1; J<=MAXXL; J++) {
          JJ=J;
          if(FACROT[J]>UEX) break;     //Oleg check this break
          TWOSR[J]=2.*J-1.+2.*SR;
          }

      if(J<=MAXXL) {
              JMX=JJ-1;
              if(JMX <=0) {
                      KK++;
                      RLEV[KK]=0;
                      MAXJ[IE]=1;
                      MAXJS[IE]=KK;
                      continue;      //Oleg check this continue
                      }
              }

       for(J=1; J<=JMX; J++) {
              AVLOG=log(FACT2)/2.0000+sqrt(ALIT*(UEX-FACROT[J]))-43.74912;

              if (AVLOG>707.) {    // this if was commented before   OLEG    // it was before 44
                     fprintf(f09,"\\par ***** E EXC.* MASS.CN EXCEEDED CCPN CAPACITY %8.3f",AVLOG);
                     s << "<p>  ***** E EXC.* MASS.CN EXCEEDED CCPN CAPACITY " << qbuf.number(AVLOG) << "</p>";
                     ALIT=0.; return;
              }

              RLEVD1= exp(AVLOG);
              RLEV[KK+J]=TWOSR[J]*RLEVD1*RLEVD1;
              }

      KK+=JMX;

      if(KK>=DIM_RLEV) {
              fprintf(f09,"\\par  ***** LEVEL DENSITY EXCEEDED %d",DIM_RLEV);
              s << "<p>   ***** LEVEL DENSITY EXCEEDED " << qbuf.number(DIM_RLEV) << "</p>";
              ALIT=0.;
              return;
      }
      s.flush();
      MAXJS[IE]=KK;
      MAXJ [IE]=JMX;
      } //continue; L6:
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FIND(double E,int &IXRSM, double *EBIN, int MEBIN, int MODE, QTextStream& s) {
//     WRITTEN BY A. GAVRON
//     SUB. TO FIND INDEX IXRSM CORRESPONDING TO NEAREST ENERGY BELOW E


      IXRSM=1;
      if(E < 0.)return;
      if((_IZR[MODE]+_INR[MODE] ) * _IZR[MODE] <=0 ) return;

      if( ! (E>=EBIN[1] && E<=EBIN[MEBIN] )) {
      _cl_KK++;

                                                //     FIND LIMITS ERROR MEANS THAT THE ENERGY RANGE OF THE LEVEL
                                                //     DENSITY TABLE IS TOO SMALL TO ACCOMODATE THE POSSIBLE DECAY
                                                //     MODES. EITHER THIS IS AN INPUT ERROR (IF YOU ARE INPUTTING
                                                //     A SPECIFIC RANGE FOR A LEVEL DENSITY ) OR SUBROUTINE RANGE
                                                //     NEEDS MODIFICATION. CONSULT ADDITIONAL COMMENT IN SUBROUTINE RANGE.
      if(_cl_KK<=10){
              fprintf(f09,"\\par  \"FIND\" limits error. E=%.2f not between %.2f ana %.2f  for  Z=%d  N=%d"
                          "\\par Consult comments in subroutine RANGE and check input double for LEVEL DENSITY",
                           E,EBIN[1],EBIN[MEBIN],_IZR[MODE],_INR[MODE]);
              s << "<p> \"FIND\" limits error. E = " << qbuf.number(E) << " not between " << qbuf.number(EBIN[1]) << " and "
        << qbuf.number(EBIN[MEBIN]) << " for Z = " << qbuf.number(_IZR[MODE]) << " N = " << qbuf.number(_INR[MODE])
              << "</p><p>Consult comments in subroutine RANGE and check input double for LEVEL DENSITY</p>";
      }

      if(_cl_KK==10){
               fprintf(f09,"\\par%s Tenth warhning. Further warnings suppressed. %s"
                           "\\par%s        Better look it.                       %s",
                           s10star,s10star,s10star,s10star);
               s << "<p>Tenth warning. Further warnings supressed.</p>";
      }


      IXRSM=MEBIN;
      if(E<EBIN[1])IXRSM=1;
      }
else  {
      IXRSM=E-EBIN[1]+1.05;
      IXRSM = min(IXRSM,MEBIN);

      while(EBIN[IXRSM]<=E) IXRSM++;

      do IXRSM--; while(EBIN[IXRSM] > E);
      }
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PACE::CHNPRB(int IPOT,int IXR, fusion_event &c, double SC,double SE,double *RLEV,double& PROBM, int *MAXJ,
        int *MAXJS, double *EBIN, int IXMIN, QTextStream& s)
{
//
//      CHNPRB(K,IXRSM[K],JC,SC,_SPART[K],
//             _RLEV[K],PROBM[K],_MAXJ[K],_MAXJS[K],_EBIN[K], AMR, IXMIN[K]);
//
//
//      IPOT - mode (n,p,a,g)
//      IXR  - from FIND
//      ISC = SPIN OF LEVEL
//      SE = PARTICLE SPIN
//      SC  = BASIC SPIN UNIT  ( 0 OR 1/2 )
//      PROBM - probability

//      AMR = A of compound
//      IXMIN = from FIND

//      SR  = MINIMUM SPIN WE CAN DECAY TO
//      STARTING FROM HIGHEST LEVEL IXRSM , WE GO LEFT (IGO=-1)
//      OR RIGHT TO LOOK FOR MAXIMUM PROBABILITY.
//
//
//     MODIFIED JULIAN SUBROUTINE
//
///   AMR -> c.Ad()
//    ISC ->  c.J
//   _ERGC -> c.Ex


static int  MAXJT[DIM_EBIN], MAXJST[DIM_EBIN];
static double  EBINT[DIM_EBIN],TL[25] ,TLSUM[25];
static int MEBINT=0;
static int IXPR=0;


long double zdum1,zdum2;

//
//      WEISSKOPF UNITS FOR GAMMA DEEXCITATION
//
const double WUE[4]={0, 6.8E-08,4.9E-14,2.3E-20};
const double WUM[4]={0, 2.1E-08,1.5E-14,6.8E-21};

extern int _cl_IZPREV;
extern int _cl_INPREV;

const double PI2=6.28318;
const double TOTHRD=0.666667;


double   SR=fabs(SC-SE);
      int IXRSM=IXR;

      int JSC=c.J;
      int JMAX=c.J;
      int MODUL=(IPOT-1)*79696;
      int AJC=c.J-1+SC;
      double ECM=0, AJR=0, SUMTL=0, DJ=0, SJ=0;
      int IGO=0, JSR=0, I, IZZ, INN, LLMAX=0, LLM=0, LMAX=0, LMIN=0, LM1, LM2;
      double /*SMIN=0, SMAX=0,*/ ALMAX=0, ALMIN=0, PRO=0;
      int MJS=0;

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
if(IPOT!=4)

   while (true) {
      ECM=c.Ex+_BE[IPOT]-EBIN[IXRSM];
      IZZ=_IZR[IPOT];
      INN=_INR[IPOT];
      TLLL(ECM,IPOT,IZZ,INN,TL,LLMAX,s);
      if(TL[1]<1.E-10) goto L58;
      TLSUM[1]=TL[1];
      LLM=LLMAX+1;
      for(I=2; I<=LLM; I++) TLSUM[I]=TLSUM[I-1]+TL[I];

      IGO =-1;                           //      JSR = J OF RESIDUAL NUCLEUS.
      JSR=JSC;
      JSR=min(JSR,MAXJ[IXRSM]);

L5:   AJR=JSR-1+SR;
      SUMTL=0.;
      DJ=fabs(AJR-AJC);
      SJ=AJR+AJC;

      ALMAX=SJ+SE;                      //      ALMAX = SUM OF SPINS
      LMAX=min(ALMAX+0.01,LLMAX);

                                 //      MINIMUM L CONTRIBUTING TO TRANSITIONS
      if(SE<DJ)         {  ALMIN=DJ-SE;  LMIN=ALMIN+0.1;}
      else if(SE>SJ)    {  ALMIN=SE-SJ;  LMIN=ALMIN+0.1;}
      else                               LMIN=0;

      if(LMIN>LMAX)goto L6;

                //  ++++++
                //  ++++++ FOLLWING IS SUBSTITUTE CODE FOR CHANNEL SPIN
                //  ++++++ OF  0  AND  1/2 .

      LM1=LMIN;
      LM2=LMAX+1;
      SUMTL=TLSUM[LM2];
      if(LM1>0)SUMTL=TLSUM[LM2]-TLSUM[LM1];
      SUMTL+=2.*SE*(SUMTL-TL[LM2]-TL[LM1+1]);


                //     SUMTL IS BEING MULTIPLIED BY  1 FOR SE=0  ALL L-S
                //                                   2 FOR SE=1  L ABOVE LMIN BELOW LMAX
                //                                   1 FOR SE=1  L=LMIN OR LMAX
                //  ++++++ END OF SUBSTITUTE CODE.
                //  ++++++
                // **********
                // ********** FOLLOWING IS ORIGINAL SUM OVER TL
                // ********** TILL -SCANNING FOR ALL L-
                // ********** IT ENABLES COUPLING OF CHANNEL SPIN
                // ********** GREATER THAN 1/2.
                // ********** FOR N,P,ALPHA,GAMMA,PRESENT CODING IS ADEQUATE.
                //  13 DO 14 L=LMIN,LMAX
                //     ABOVE LOOP OFTEN STARTS AT ZERO. IF THIS BOTHERS YOUR
                //     COMPILER, SHIFT THE INDICES BY ONE UNIT UP.
                //     DL=ABS(AJC-L)
                //     SL=AJC+L
                //     S1=AMAX1(SMIN,DL)
                //     S2=AMIN1(SMAX,SL)
                //     IDS=S2-S1+1.01
                //     SUMTL=SUMTL+IDS*TL(L+1)
                //  14 CONTINUE
                //      SCANNING FOR ALL L IN INTERVAL
                // ********** END GENERALIZED CHANNEL SPIN COUPLING CODE.
                // **********
                //      ORDINAL NO. OF THIS ENERGY LEVEL
                //PAR  PRO=RLEV(MJS)*SUMTL*.5
                //     IGNORE PARITY FACTOR.DOES NOT EFFECT RESULTS

      MJS=max(1,MAXJS[IXRSM]+JSR-MAXJ[IXRSM]);   // Oleg 06/13/2011
//      MJS=MAXJS[IXRSM]+JSR-MAXJ[IXRSM]
      if (SUMTL<1E-20) PRO=0.;
      else  {
            if( (log(RLEV[MJS])+log(fabs(SUMTL))) > -87)  PRO=RLEV[MJS]*SUMTL;
            else  PRO=0.;
            }


      if (_DISCPR <= 0. || PROBM <= 0.) zdum1=0.;
      else {
           zdum1=logl(_DISCPR);
           zdum2=logl(PROBM);
           if ((zdum1+zdum2)<-87.)  zdum1=1.E-38;
           }

      if(PRO<=zdum1)goto L6;
//      SEE REMARK ABOUT DISCPR IN MAIN PROG.
      _IPROB++;
                        //      IPROB IS TOTAL NO. OF PROBABILITIES CALCULATED
                        //     IF(IPROB.GT.5555)GO TO 60
      _PROB[_IPROB]=PRO;
      _MJ[_IPROB]=MJS+MODUL;
      if(PRO>=PROBM) {
              PROBM=PRO;
              JMAX=JSR;
              }             //      SEARCH IN ALL EX AND J DIRECTIONS FOR PROB.GE..001PROBMAX


L11:  if(JSR==1)goto L10;
      JSR+=IGO;

L9:   if(JSR>MAXJ[IXRSM])goto L58;
      goto L5;




L6:   if(JSR==JSC)goto L11;

      if(IGO==1)goto L58;

L10:  IGO=1;
      JSR=JSC+1;
      goto L9;

L58:  if(IXRSM==IXMIN)return;

      IXRSM--;
      JSC=JMAX;


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WW    POT = 4    WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
} else {
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW



while(true) {
      double AMRT=pow(c.Ad(),TOTHRD);
      double AMRTT=pow(AMRT,2);
      double PAMRT=PI2*AMRT;

      _GE1=_FGE1*PAMRT*WUE[1];
      _GM1=_FGM1*WUM[1]*PI2;
      _GE2=_FGE2*PI2*AMRTT*WUE[2];
      _GM2=_FGM2*PAMRT*WUM[2];

      ECM=c.Ex-EBIN[IXRSM];
      double ECUBE=ECM*ECM*ECM;
      double EFIVE=ECUBE*ECM*ECM;
      double GGE1=_GE1*ECUBE;
      double GGM1=_GM1*ECUBE;
      double GGE2=_GE2*EFIVE;
      double GGM2=_GM2*EFIVE;
//PAR  PROBA=(GGE1+GGM1)*.5
//PAR  PROBB=(GGE2+GGM2)*.5
      double PROBA=GGE1+GGM1;
      double PROBB=GGE2+GGM2;
//     IGNORE PARITY
      double PROBC=PROBA+PROBB;
      int L=JSC-3;
      int JR;
      double PROB1=0;

for(I=1; I<=5; I++) {
      JR=I+L;
      AJR=JR-1+SR;

      if(JR<=0 || JR>MAXJ[IXRSM] || (AJC<0.1 && I==3)  ) continue;

             if    (I==1||I==5)                                 PROB1=PROBB;
      else   if( ( (I==2||I==4) && AJC<0.1 && AJR<0.1 )   ||
                   (      I==3  && AJC<0.6 ))                   PROB1=PROBA;
             else                                               PROB1=PROBC;

      MJS=MAXJS[IXRSM]+JR-MAXJ[IXRSM];

      PRO=RLEV[MJS]*PROB1;
      if (_DISCPR<=0.||PROBM<=0.) zdum1=0.;
      else {
           zdum1=logl(_DISCPR);
           zdum2=logl(PROBM);
           if ((zdum1+zdum2)<-87.) zdum1=1.E-38;
           }

      if(PRO <= zdum1)continue;

      _IPROB++;

      if(_IPROB>5555) {
             fprintf(f09,"\\par DIMENSION OF PROB NOT HIGH ENOUGH - ****************");
             s << "<p> DIMENSION OF PROB NOT HIGH ENOUGH - ****************</p>";
             return;
             }

      _MJ[_IPROB]=MJS+MODUL;
      _PROB[_IPROB]=PRO;
      PROBM=max(PROBM,PRO);
     }

  if(IXRSM==IXMIN) break;
  IXRSM=IXRSM-1;
  JSC=JMAX;
  }

      double BARR=_FBARR[c.J];
      _LFISS=0;
      if(BARR>_BARMAX)return;

                //********** ATTENTION *******
                //     NO FISSION CALCULATED IF BARRIER HIGHER THAN 30 MEV
                //     AT THIS ANGULAR MOMENTUM

      if(c.Ex<BARR+EBIN[1])return;


      if( _cl_IZPREV!=_IZR[4] || _cl_INPREV!=_INR[4] ) {
              _cl_IZPREV=_IZR[4];  _cl_INPREV=_INR[4];
              LEVDEN(_IZR[4],_INR[4],1,c.Ex+2.5, MAXJT, MAXJST, EBINT, MEBINT, RLEVT,IXPR,_BE[4],_ARATIO,s);
              }

      FIND(c.Ex - BARR,IXRSM,EBINT,MEBINT,1,s);
      double SUM=0, R=0;
      int N, NLEV;

   for(I=1; I<=IXRSM; I++) {
      N=IXRSM-I+1;
      if(c.J>MAXJT[N])continue;
      NLEV=MAXJST[N]-MAXJT[N]+c.J;
      R=RLEVT[NLEV];
      SUM+=R;
      if(R<0.001*SUM)break;
      }

       PRO=SUM;
                              //     IGNORE PARITY
      _LFISS=1;
      _IPROB++;
      if(_IPROB>5555){
            fprintf(f09,"\\par DIMENSION OF PROB NOT HIGH ENOUGH - ****************");
            s << "<p> DIMENSION OF PROB NOT HIGH ENOUGH - ****************</p>";
             s.flush();
            return;
            }
      s.flush();
      _MJ[_IPROB]=_IMODL;
      _PROB[_IPROB]=PRO;
      return;
     }
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void SEARCH(fusion_event *a, fusion_event &c, int &MODE, QTextStream& s)
//void SEARCH(int &IZ, int &IN, double &ERG, int &J, int  IZC, int INC, int &MODE)
 {
//
//     ORIGINAL JULIAN SUBROUTINE
//
//     ***********
//
//     SUBROUTINE  DETERMINES FINAL RESIDUAL NUCLEUS
//      SEARCH(_IZCS[I],_INCS[I],_ERGCS[I],_JCS[I],IZC,INC,MODE,IRAND)

// ATTENTION with MAX*[][]  !  they was transposed!

int I, JF, IS, MJI=0, IBIN=0, MEB=0;

double R=RANF();

                //      THE SEARCH ALGORITHM IS TO KEEP HALVING THE VECTOR AND
                //      DECREASING LIMITS TILL INDEX IS FOUND. THE INDEX OF THE
                //      DETERMINES THE DECAY MODE

if(R > _PROB[_IPROB])  I=_IPROB;
else {
      if(R <= _PROB[1]) I=1;
      else {
              IS=0;
              JF=_IPROB;

        L2:   JF=JF/2+1;
              if(JF==2)JF=1;
              IS = min(IS+JF,_IPROB);

        L3:   if(R>_PROB[IS])goto L2;

        L1:   JF=JF/2+1;
              if(JF==2)JF=1;
              IS=max(IS-JF,1);
              if(JF!=1)goto L3;

              if(R<=_PROB[IS])goto L1;
              I=IS+1;
           }
      }
      MODE=_MJ[I]/79696+1;

      if(MODE>=5) {
                a->Ex=-1.; // it means FISSION
                return;
                }

//      IBIN IS PLACE IN 79696 LONG VECTOR OF ONE SPECIFIC MODE.
      a->init(c.Z-_IZPART[MODE],c.N-_INPART[MODE]);

      MJI=_MJ[I];
      IBIN=((MJI)%(79696));
      if(IBIN==0)IBIN=79696;
      MEB=_MEBIN[MODE];
      for(I=1; I<=MEB; I++) if(_MAXJS[MODE][I] >= IBIN)goto L10;

      fprintf(f09,"\\par SEARCH ERROR. MODE=%d MEBIN=%d MAXJS=%d IBIN=%d",MODE,MEB,_MAXJS[MODE][MEB],IBIN);
      s << "<p>  SEARCH ERROR. MODE = " << qbuf.number(MODE) << " MEBIN = " << qbuf.number(MEB) << "MAXJS = " <<
           qbuf.number(_MAXJS[MODE][MEB]) << "IBIN = " << qbuf.number(IBIN) << "</p>";
      s.flush();
      return;

L10:  a->J  = _MAXJ[MODE][I]-_MAXJS[MODE][I]+IBIN;
      a->Ex = _EBIN[MODE][I];

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void MJRAN(fusion_event &csi, fusion_event &csf, int MODE, double EP, QTextStream& s)

{

//void MJRAN(int JI,int JF,int MI,int &MF,int NP,double EP,int IZ,int IN)
//     MJRAN(csi.J,csf.J,csi.MJ,csf.MJ,MODE,EP,csi.Z,csi.N);
//void MJRAN(fusion_event &csi, fusion_event &csf, int MODE, double EP);

//     WRITTEN BY A. GAVRON
//     THIS SUBROUTINE DETERMINES THE M-STATE DISTRIBUTION OF THE
//     ANGULAR MOMENTUM I.E., THE ORIENTATION IN SPACE OF THE COMPOUND
//     NUCLEUS.
//     ********* WARNING -  THIS IS DONE IGNORING THE SPIN OF THE
//     EMITTED PARICLE SINCE IT IS NOT JUDGED TO HAVE A SIGNIFICANT
//     EFFECT ON THE RESULTS.
//     THE REASON FOR INCLUDING M-STATE DISTRIBUTION IS TO OBTAIN
//     CORRECT ANGULAR DISTRIBUTIONS FOR PARTICLES AND RESIDUAL NUCLEI.
//
//     THE INITIAL M-STATE DISTIBUTION IS ASSUMED TO BE M=0 FOR ALL
//     EVENTS PRODUCED BY A COMPOUND NUCLEUS REACTION.
//     FOR A FRAGMENT PRODUCED IN A DEEP INELASTIC COLLISION, M=J
//     IS ASSUMED; THE Z AXIS IS THEN TAKEN PERPENDICULAR TO THE
//     REACTION PLANE AND THE DIRECTION OF MOTION OF THE FRAGMENT IS
//     THE X AXIS.
//     MJCS(I) IN THE MAIN PROGRAM IS ACTUAL M VALUE - NOT AN INDEX
//             AS USED FOR JCS(I)  (WHICH IS J+1 OR J+0.5)
//             AGAIN NOTE THAT WE USE INTEGER M'S AND NEGLECT SPIN.

      double TL[25],TM[50];

//      common/XQANG/SUM[300][10],MDIR,CT;

 //     SUM IS TABLE OF INTEGRALS FOR CHOSING ANGLE FOR A GIVEN SPHERICAL
//     HARMONIC
      int LMAX=0;
      TLLL(EP,MODE,csi.Z,csi.N,TL,LMAX, s);
//     FROM TRANSMISSION COEFF VECTOR WE SELECT THE L VALUE THAT WAS
//     RESPONSIBLE FOR THE TRANSITION FROM  JI  TO  JF.
// - LMAX IS ACTUAL L  (NOT  L+1)

      int I;

      double X, FAC;

      int L1=abs(csi.J-csf.J);
      int L2=csi.J+csf.J-2;
      L2=min(L2,LMAX);
      L1++;
      L2++;
      int LA=L1+1;
      int LL=L1;
  if(LA<L2){//----------------------------------
      for(I=LA; I<=L2; I++) TL[I]+=TL[I-1];
      double FAC=1./TL[L2];
      double X=RANF();

      for(I=L1; I<=L2; I++) {
              LL=I;
              if(X<=TL[I]*FAC)break;
               }
 }//----------------------------------
      int L=LL-1;            //     M-S ARE ACTUAL VALUES AND NOT INDICES ( I.E. M+1)
                          //     L IS NOW KNOWN- WE PROCEED TO CHOOSE  M.
      int MA=csf.J-1+csi.MJ;
      MA=min(MA,L);
      MA=-MA;

      int MB=csf.J-1-csi.MJ;
      MB=min(MB,L);
      int NM=MB-MA+1;

      double FJI=csi.J-1;
      double FJF=csf.J-1;
      double FL=L;
      double FMI=csi.MJ;
      int M=1;

 if(NM>1) {//----------------------------------

      for(I=1; I<=NM; I++) {
              M=MA+I-1;
              double FM=M;
              double FMF=-FMI-FM;
              TM[I]=pow(C3J(FJI,FL,FJF,FMI,FM,FMF),2);
              }

              //     M IS CHOSEN ACCORDING TO WEIGHT OF  3-J  (C.G.) SYMBOL SQUARED
              //L997: format(1X,6F8.3);

      for(I=2; I<=NM; I++) TM[I]+=TM[I-1];

      FAC=1./TM[NM];
      X=RANF();

      for(I=1; I<=NM; I++) {
              M=I;
              if(X<TM[I]*FAC)break;
              }

 } //----------------------------------

      M+=MA-1;
      csf.MJ=csi.MJ+M;
      M=abs(M);
      int N=LL*(LL-1)/2+M+1;
                                //     FIND SUM TABLE FOR THIS L-M COMBINATION
      X=RANF();
      for(int ICT=1; ICT<=10; ICT++) {
              _CT=ICT;
              if(X<=_SUM[ICT][N])break;   // already transposed
              }                          //      CHOOSE ANGLE COSINE
      _CT=(_CT-RANF())*0.1;
      X=RANF();
      if(X>0.5)_CT=-_CT;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double C3J(double FJ1,double FJ2,double FJ3,double FM1,double FM2,double FM3)
{
//     OBTAINED FROM AVRI FRAENKEL.
//      common/FAC/FA[203];

double *FACT;
double C3J=0;

FACT=&FA[1];   //EQUIVALENCE [FA[2]][FACT[1]];

int M1, M2, M3, M4, M5, M6, M7, K, NJ, IZ, II;
double B2, B3, B4, B5, B6;
double A3 = 0.;
double G1 = FJ1+FJ2;
double G2 = FJ1-FJ2;
double G3 = FJ3-G2;
double G4 = G1-FJ3+.01;
double G5 = G2+FJ3;
double G6 = G1+FJ3;
double A = FJ3-fabs(G2)+.01;


if(A<0||G4<0.)return 0;


    M1 = G4;
    M2 = G5+.01;
    M3 = G3+.01;
    M4 = G6+1.01;

double A1 = FACT[M1]+FACT[M2]+FACT[M3]-FACT[M4];
double H1 = FJ1+FM1;
double H2 = FJ2+FM2+.01;
double H3 = FJ3+FM3;
double P1 = FJ1-FM1+.01;
double P2 = FJ2-FM2;
double P3 = FJ3-FM3;
double P4 = G5-P1+.02;
double P5 = G3-H2+.02;
      M1 = H1+.01;
      M2 = H2;
      M3 = H3+.01;
      M4 = P1;
      M5 = P2+.01;
      M6 = P3+.01;
      M7 = M1-M5;
double A2 = (A1+FACT[M1]+FACT[M2]+FACT[M3]+FACT[M4]+FACT[M5]+FACT[M6])/2.;
double SIGN=-1.;
      NJ = G4+1.;


for(II=1; II<=NJ; II++) {
      IZ=II-1;
      SIGN=-SIGN;
      B2 = G4-IZ;
      B3 = P1-IZ;
      B4 = P4+IZ;
      B5 = P5+IZ;
      B6 = H2-IZ;
      if(B3<0 || B4<0 || B5<0 || B6<0) continue;
      M1=IZ;
      M2=B2;
      M3=B3;
      M4=B4;
      M5=B5;
      M6=B6;
      A3+=SIGN/(exp(FACT[M1]+FACT[M2]+FACT[M3]+FACT[M4]+FACT[M5]+FACT[M6]-A2));
      }

      K=M7/2;

      if(M7 - 2*K == 0) C3J= A3;
      else              C3J=-A3;

      return C3J;
      }



