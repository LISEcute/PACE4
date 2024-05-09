#include "include/Constant.h"
#include <time.h>
#include <stdio.h>
//#include <dos.h>
#include <math.h>
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define dr13 (1./3.)
#define pow_dr13(A)    pow(A,dr13)
#define pow_mdr13(A)   pow(A,-dr13)

//#include <vcl.h>

double GETSEED(void);  //  return ISEED
double CLOCK(); //return I
void FACLOG() ;
void PLM() ;
void COMPOS(int &IZ,int &IN,double &ENERGY,double &VZC,int &MAXC,
                double &SIGMA,double &ECLOSS,int INPUTp);
extern double MASSES(int IZ, int IN);
double BASS (int IA1,int IZ1,int IA2,int IZ2, double E);
bool BARFIT (int IZ,int IA,int IL, double &BFIS, double &SEGS, double &SELMAX);
void LPOLY (double X, int N, double *PL);       // base = 1
double FISROT(double A, double Z, double AL, double BARFAC);
double YRAST(double A, double Z, double AL);
void OUTEM (int ICTL, int MODE, double EEML, double AEML);
void OPTPOT (int I, int IN, int IZ, double &RR, double &AR, double &RI, double &AI);
void TLLL (double &ER, int NP, int IZ, int IN, double *TL, int &LLMAX);
void COMPND (int IZ, int IN, double ENERGY, int &MAXC);
int ISUM (int N[17][33][28],int I,int J,int L1,int L2);
char* LMNT (int IZ);
void OUTRES (int ICTRL, int LP, int NP);

extern double pow_int(double par, int power);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
//  GLOBAL from DIALOG
int NCASC, IDIST, INPUTp, FISSBR, NOSHL;
double FYRST,  FACLA;
int IZP,IAP,IZT,IAT,LMINN;
double SP,ST,QCN,ELAB,EXPSIG,AGRAZ,ELOSS;



//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww

//------COMMON /EM/
int  NT[5], N[4][63][19], NG[63], NW[4][5][19], NGW[5],
     NGB[63], NGBW[5], NGT[63], NGTW[63];
double  EW[4], EW1[5], EW2[4];
//----------------------

//-------COMMON /XQANG/
double SUM[301][11];
int MDIR;
double CT;  // COSTH
//---------------------

//-------COMMON /MASS/
double FLMASS,  BARFAC, ARATIO, FBARR[111];
int IZPART[5]={0,   0,1,2,0};
int INPART[5]={0,   1,0,2,0};
int IZR[5], INR[5], NFILE, LFISS;
//---------------------


long double FA[204];                // COMMON /FAC/
double SIGJC[101], AJCN[151],ALEV[101];

//------- common /pshel/
int inshel,iashel,izshel;
double pshel;
//

double  DSIG[19], SIGMA;
double clshft[4];
double EMAX[5], EMIN[5], BE[5];
double PROB[2223],  GE1, GE2;
int IPROB, MJ[2223], IMODL, MAXC;
double EROT[111];


//------COMMON /SCH/-------------------
int MAXJ[287][5],  MAXJS[287][5], MEBIN[5];
double EBIN[287][5], RMASS[5];
//-------------

//------COMMON /OUT/ -------------------
int IZCS[9997], INCS[9997], JCS[9997], MJCS[19997];
double ERGCS[9997], RLEV[19997][5];
//-------------

//------COMMON /GAM/
double FGE1, FGM1, FGE2, FGM2, ERGC, DISCPR,ERGMIN, EDL[37][7];
int IADEF, NDISC/* =ND*/, NDL[7], IZDL[7], INDL[7], JDL[37][7];
//---------------

//--------COMMON /XQTT/
double  UX, EX, TT;
//---------------

//--------COMMON /RANDOMM/
int IRX;
//---------------

//------- COMMON /DAT/
double SPART[5]={0.,0.5,0.5,0.0,0.0};
double PMASS[5]={0.,8.07169,7.28922,2.42494,0.0};
//---------------

//------- COMMON /XQRES/
int IZ, IN, ITRAC;         // IR, IQ, ENERGY, VZC, ITRAC, IDIST
double ENERGY, VZC;        // EX, VZ
//---------------

//------  COMMON /POT/
double ATT[6], AV[5], AW[5], AWV[5], ASO[5], AW1[5], EP[37],
        ECUT, ETHRS, RIV, AIV, RSOT, ASOT,  RCX, VQS[4], VQ0[4];
int  IRAD, IMT,  IMESH;
//-------------------


//-------- COMMON /Q10/
double PC, PM, TC, TM, V0, W0, /*RCX,*/ WVOLF, AT[6], RMATCH;
int IIM, /*IMESH,*/ JUC;
//-------------------

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
FILE *f10;
FILE *f04;
const char *f10_name="result.out";
const char *f04_name="result.opt";
char *s35x="         ";
char *s100m="---------------------------------------------------------------------\n";

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWww
void start(void) {

//#include "start.h"
int NSPC[5][52]; //4,51)
char *OPTFIL="results.opt";
int ISEED, IRX,TT1, JTIME1;
double TL[31];

time_t t;
/*              CHARACTER*8 FILNAMS
                CHARACTER*12 MASSFIL,LISFIL,INFIL,EVFIL,OPTFIL,PARFIL
                CHARACTER*4 IEXT(5)
//         }//end VAX I/O
                EXTERNAL XDATA
                DIMENSION PROBM(4), EXRSM(4), IXRSM(4), IXPR(4), IXMIN(4), TL(30),, IBUF2(9), IMGZ(9996)
//          VAX I/O
                DATA IEXT/'.INP','.OUT','.OPT','.EVT','.PAR'/
//
                REAL*8 PROB,SPROB,TPROB
//              CALL ERRSET(219,1,1,0,0,220)
//          VAX I/O
                READ 9002,FILNAMS
9002		    FORMAT(A8)
9003            FORMAT(A8,A4)
                WRITE(LISFIL,9003)FILNAMS,IEXT(2)
                WRITE(INFIL,9003)FILNAMS,IEXT(1)
                WRITE(EVFIL,9003)FILNAMS,IEXT[4]
                WRITE(OPTFIL,9003)FILNAMS,IEXT[3]
                WRITE(PARFIL,9003)FILNAMS,IEXT[5]
                OPEN(UNIT=9,ACCESS='SEQUENTIAL',STATUS='OLD',ERR=8009,  FILE=INFIL)
                OPEN(UNIT=2,ACCESS='SEQUENTIAL',STATUS='SCRATCH',ERR=8002,  FORM='UNFORMATTED',FILE=EVFIL)
                OPEN(UNIT=3,ACCESS='SEQUENTIAL',STATUS='OLD',ERR=8003,  READONLY,FILE='[signorini.pace2_pu]MASSES.DAT',SHARED)
                OPEN(UNIT=4,ACCESS='SEQUENTIAL',STATUS='NEW',FILE=OPTFIL)
//           CC
*/
f10=fopen(f10_name,"wt");
f04=fopen(f04_name,"wt");

int IDIM=19696;
int ISTOP=1;
//          IBM CALL JSTIME(JTIME1)
//          TMU CALL GETTR (TIME1)
ISEED=GETSEED();
IRX=ISEED*2+1;
TT1=time(0);
//                CALL LIB$INIT_TIMER
JTIME1 = CLOCK();
//              CALL RNDMZ
//OLEG   temporary-->   PLM();

int MODES=4;
int IMODL=4*19696;
BE[4]=0;
double ECLOSS=0;
EMIN[1]=0.01;
EMIN[4]=0.;
//                WRITE(10, 850 )
int IBOUND=0;
int INOTGS=0;
int INOTG1=0;
double SNOTGS=0;
double ENOTGS=0;
int NODES=MODES-1;
//
//READ(9,940)NCASC,IDIST,INPUTp,FYRST,FISSBR,ARATIO,MDIR,FACLA, ITRAC,NOSHL,
//?????????NTIME,clshft(2),clshft(3),pshel,iashel,izshel
iashel=0;
clshft[1]=0;
clshft[2]=0;
clshft[3]=0;
inshel=iashel-izshel;
//
int IDUMP=0;         if(INPUTp<0 || ITRAC<0) IDUMP=1;

INPUTp=abs(INPUTp);
ITRAC =abs(ITRAC);
NCASC = min(9996,NCASC);

//              ARATIO IS (LITTLE-A-SADDLE)/(LITTLE-A-G.S.)
//              DEFAULT FOR NO INPUTp IS 1.000

if(ARATIO==0) ARATIO=1;
double ACASC=(ACASC < 100 ? 100 : NCASC);


//              DISCPR IS DISCRIMINATOR ON PROBABILITY VECTOR (PROB) TO PREVENT
//              SEARCH FOR HIGHLY UNLIKELY EVENTS. CHANGE TO 0.001 IF YOU
//              WANT RARE EVENTS AND ARE USING GOOD STATISTICS.

double DISCPR=1./ACASC;


//              COMPOS FORMS INITIAL SPIN DISTRIBUTION DEPENDING ON INPUTp
COMPOS(IZ,IN,ENERGY,VZC,MAXC,SIGMA,ECLOSS,INPUTp);


int IADEF=2;
if(IN>=77&&IN<=110) IADEF=1;
if(IN>138) IADEF=1;
if(IADEF==2&&FACLA<0.)fprintf(f10,"\n%sSpherical nucleus level density ******",s35x);
if(IADEF==1&&FACLA<0.)fprintf(f10,"\n%s Deformed nucleus level density ******",s35x);

//              IADEF=2 IS ASSUMED SPHERICAL NUCLEI FOR LEVEL DENSITIES.
//              IADEF=1 IS CONSIDERED A DEFORMED NUCLEUS.
//              THE DECISION IS MADE ONLY FOR THE Z,N OF THE INITIAL NUCLEUS.
int IZMAX=IZ;
int INMAX=IN;
double FLMASS=MASSES(IZ,IN);
double BARFAC=1.;
int IA=IZMAX+INMAX;
double A=IA;
IZ=IZMAX;
double SELMAX, SEGS,BRR;
BARFIT(IZ,IA,0,BRR,SEGS,SELMAX);

if(FISSBR!=0.) {
     BARFAC=FISSBR/BRR;
     if(FISSBR<=0.){
         BARFAC=-FISSBR;    //  NEGATIVE FISSBR MEANS IT IS THE FACTOR IN ABS. VALUE
         FISSBR=BRR*BARFAC; //  BARFAC MULTIPLIES THE SIERK BARRIER TO GET YOUR OWN.
         }
fprintf(f10,"%sInput fission barier                         %6.2f  MeV\n", s35x,FISSBR);
fprintf(f10,"%sSierk barrier will be multiplied by          %6.2f\n", s35x,BARFAC);
      }
else
fprintf(f10,"%sZero spin fission barrier                    %6.2f\n", s35x,BRR);

fprintf(f10,"%slittle-a sub-f / little-a sub-gamma          %6.2f\n", s35x,ARATIO);


if(FACLA <0)
fprintf(f10,"%sGilbert and Cameron little-a selected         *****\n", s35x);

if(FACLA==0){
fprintf(f10,"%sDefault value of FACLA substituted            *****\n", s35x);
             FACLA=7.5;
             }

if(FACLA>0)
fprintf(f10,"%slittle-a mass                                  %4.1f\n", s35x,FACLA);

if(FYRST <0)
fprintf(f10,"%sGilbert and Cameron yrast line used           *****\n", s35x);

if(FYRST==0) {
fprintf(f10,"%sFYRST = 0.   Changed to default value         *****\n", s35x);
FYRST=1;
}

if(FYRST >0) {
fprintf(f10,"%sSierk diffuse surface nucleus yrast line used\n", s35x);
fprintf(f10,"%sMultiplied by factor of                       %5.2f\n", s35x,FYRST);
}

if(NOSHL==1)
fprintf(f10,"%sLysekil masses with no shell correction used  *****\n", s35x);

double TEMP=sqrt(ENERGY*10./A);
double A3=pow_dr13(A);
double RD=A3+1.;
double VCOULP=IZMAX/RD;
RD=A3+1.587;
double VCOULA=2.*IZMAX/RD;
EMAX[1]=16.;
EMAX[4]=16.;
//              ***************************
EMAX[3]=VCOULA*1.1+5.*TEMP;
EMIN[2]=0.4*VCOULP;
EMIN[3]=0.5*VCOULA;
EMAX[2]=VCOULP*1.1+5.*TEMP;

//          THE EMIN ARE THE SMALLEST ENERGIES FOR WHICH RELIABLE CROSS SECTIONS
//              ARE CALCULABLE.  THE FORMULAE WERE OBTAINED FROM EXPERIENCE.
//              EMAX ARE MAXIMUM EMISSION ENERGIES SO AS NOT TO WASTE TIME
//              COMPUTING UNLIKELY EVENTS
for(int IKQ=1; IKQ<=4; IKQ++)
   for(int JKQ=1; JKQ<=51; JKQ++)
               NSPC[IKQ][JKQ]=0;

fprintf(f10,"\n    Ecm range for   neutron   proton    alpha    gamma   (MeV)\n");
fprintf(f10,"                   %8.2f %8.2f %8.2f %8.2f    min\n",EMIN[1],EMIN[2],EMIN[3],EMIN[4]);
fprintf(f10,"                   %8.2f %8.2f %8.2f %8.2f    max\n",EMAX[1],EMAX[2],EMAX[3],EMAX[4]);

fprintf(f10,"\n Number of cascades %i\n Internal probabi1lity discriminator of program set to %8.4f\n",NCASC,DISCPR);

int I=4;        // CHECK!!!

if(IDIST!=0) OUTEM (1,I,A,A);   //   OUTEM IS PARTICLE EMISSION OUTPUT SUBROUTINE

fprintf(f10,"Optical model parameters for light emitted particles outputted to \"%s\"\n",OPTFIL);
fprintf(f10,"Coulomb Energy Shift for Transmition Coefficient for proton %6.2f"
            " and for alpha %6.2f\n",clshft[2],clshft[3]);

if(iashel!=0)
fprintf(f10,"Shell effect was included in A=%i Z=%i and the correction for an=%.3f\n",iashel,izshel,pshel);

double RR,AR,RI, AI,EDUM;
int IDUM;
//***************************
for(I=1; I<=3; I++) {
               // READ (9,820) (AV(J),J=1,3),(AW(J),J=1,3),(AWV(J),J=1,3),RR,AR,RI,AI,RIV,AIV,RCX
                if(AV[1]==0 && AW[1]==0) OPTPOT (I,IN,IZ,RR,AR,RI,AI);
                IMESH=100;;;;
                IRAD=1;
                if(RR>1.8 || RI>1.8 || RIV>1.8) IRAD=0;
                ECUT=0;
                ETHRS=0;
                RSOT=0;
                ASOT=99;
                if(AR==0) AR=99;
                if(AI==0) AI=99;
                if(AIV==0) AIV=99;
//              ZERO DIFFUSENESS IS DEFAULT FOR NO POTENTIAL, BUT PROGRAM TRIES TO
//              CALCULATE ANYWAY. SO A-S HAVE FINITE VALUE TO PREVENT BLOWUP
                IMT=1;
//              IMT=1  IS SAXON WOOD SURFACE DERIVATIVE
//              IMT=2 (CHANGE PROGRAM) IS GAUSSIAN SURFACE POTENTIAL
                ATT[1]=RR;
                ATT[2]=AR;
                ATT[3]=RI;
                ATT[4]=AI;
               for(int J=1; J<=4; J++) { AW1[J]=ASO[J]=0;}
                TLLL(EDUM,IDUM,IZMAX,INMAX,TL,IDUM);
               }


                if(NCASC==0) return;


//  move to end
fclose(f10);
fclose(f04);
}//////////////////////////////////   END
/*

//              TLLL PREPARES ALL N, P, AND ALPHA TRANSMISSION COEFFICIENTS
//              IN ADVANCE.
//              THE CHANGE IN TRANSMISSION COEFFICIENTS FOR SECOND AND SUB-
//              -SEQUENT DECAYSIS SIMULATED BY SHIFTING THE ENERGY AT WHICH
//              THE TL-S ARE CALCULATED ACCORDING TO THE SHIFT OF THE COULOMB
//              BARRIER OF THE SECONDARY RESIDUAL NUCLEI.
//              TL-S AT THIS ENERGY ARE CALCULATED BY INTERPOLATION
//
//              **************************
                READ (9,980) FGE1,FGM1,FGE2,FGM2
if(FGE1==0.)  {
      if(A<=100.){FGE1=.00008;    FGE2=4.8;  FGM1=.025;   FGM2=.0195;}
 else if(A<=126.){FGE1=.000014;   FGE2=5.9;  FGM1=.01;    FGM2=.00088;}
 else if(A<=150.){FGE1=.000046;   FGE2=7.7;  FGM1=.007;   FGM2=.058;  }
 else            {FGE1=.000011;   FGE2=9.0;  FGM1=.010;   FGM2=1.2;   }
     }

             100 WRITE(10, 970) FGE1,FGM1,FGE2,FGM2
                if(INPUTp==2)goto 101
                WRITE(10, 830)
                WRITE(10, 840) (AJCN[I],ALEV[I],I=1,MAXC)
             101NDISC=1
//              DISCRETE LEVELS (IF DESIRED) FOR VARIOUS LEVEL DENSITIES
//              ARE NOW READ IN. THE NUMBER OF FINAL NUCLEI
//              THAT HAVE DISCRETE LEVELS IS NDISC.
//              NOTE *****  NDISC SHOULD NOT BE MORE THAN 6.
//              *****  THE MAXIMUM NUMBER OF LEVELS FOR EACH NUCLEUS IS 36
             110READ (9,990,END=120) NDL(NDISC),IZDL(NDISC),INDL(NDISC)
                NLL=NDL(NDISC)
                if(NLL>0&&NLL<=36) goto 130
             120NDISC=NDISC-1
                goto 140
             130READ(9,1000)(EDL(I,NDISC),JDL(I,NDISC),I=1,NLL)
                WRITE(10,1010)IZDL(NDISC),INDL(NDISC),
                1       (EDL(I,NDISC),JDL(I,NDISC),I=1
               1,NLL)
                NDISC=NDISC+1
                goto 110
             140IF (INPUTp==2) goto 150
//              INPUTp = 1  - PROJECTILE + TARGET
//               =2  - COMP. NUC. , EX AND J FIXED
//              = 3  - COMP. NUC. , EX FIXED, J DISTRIBUTED
//              = 4  - COMP. NUC. ,J DISTIBUTION CALCULATED INTERNALLY
             150ID=1
//              PREPARE TABLES FOR RECOILING NUCLEUS.
                if(IDIST!=0) CALL MOMENT (0,0.,0.,0,VZC,NCASC,ID)
//
                if(INPUTp==2) goto 210
//           TH E SPIN DISTRIBUTION OF THE STARTING EXCITED NUCLIDE HAS BEEN
//              DETERMINED.  THE DETERMINATION OF THE PROBABILITIES OF
//              THE PATHS OF DEEXCITATION WILL BE DETERMINED IN THE FOLLOWING
//              WAY
//              1) USING RANDOM NUMBERS, NCASC CHANNELS ARE ESTABLISHED FOR
//              THE COMPOUND NUCLEUS.
//              2) THE HIGHEST A IN THE ARRAY IS DETERMINED.
//              3) THE FIRST NUCLEUS OF THIS A IN THE ARRAY IS DETERMINED.
//              4) THE HIGHEST ENERGY FOR THIS NUCLEUS IS DETERMINED.
//              5) THE FIRST SPIN OF THIS NUCLEUS AND ENERGY IS FOUND.
//              6) THE PROBABILITY MATRIX FOR THIS NUCLEUS ENERGY AND SPIN
//              IS DETERMINED.
//              7) BY MEANS OF RANDOM NUMBERS, THE DEEXCITATION PATH IS CHOSEN
//              FOR ALL MEMBERS OF THE ARRAY OF THAT NUCLIDE, ENERGY AND
//              SPIN.  THE NEW Z, N, ENERGY, AND SPIN ARE SUBSTITUTED FOR
//              THE OLD VALUES IN THE ARRAY.
//              8) THE NEXT SPIN FOR THE (OLD) NUCLIDE AND ENERGY IS LOCATED.
//              IF IT EXISTS, goto STEP 7.
//              9) THE NEXT HIGHER ENERGY FOR THAT (OLD) NUCLIDE IS LOCATED.
//              IF IT EXISTS, GOT TO STEP 5.
//              10) THE NEXT NUCLIDE FOR THAT (OLD) A IS LOCATED.
//              IF IT EXISTS, goto STEP 4.
//              11) goto STEP 2.
//
//           SPROB   IS UNNORMALIZED SUM OF PROBABILITIES
//           PROB    (EQUIVALENT - PROB) IS UNNORMALIZED PROBABILITIES
//           MODES   IS NUMBER OF DIFFERENT PARTICLES (AND GAMMAS)
//              CONSIDERED FOR EVAPORATION
//           IZPART  AND  INPART  ARE Z AND N OF EMITTED PARTICLES
//              THE VECTOR ASSOCIATED WITH EACH EVENT IS AS FOLLOWS
//
//              ERGCS[J]- THE EXCITATION ENERGY
//              IZCS[J]- THE ATOMIC NUMBER ( Z )
//              INCS[J] - THE NUMBER OF NEUTRONS ( N )
//              JCS[J]  - THE SPIN INDEX  I.E.  1,2,3,   ETC.
//              SPIN = SPIN INDEX - 1  FOR EVEN A
//              SPIN = SPIN INDEX - 1/2 FOR ODD A
//              MJCS[J] - THE SPIN PROJECTION  =  0,1,2,   SPIN.
//              WE NEGLEGT 1/2 INTEGER VALUES FOR M *******
//
//           SPIN DISTRIBUTION OF INITIAL NUCLIDE NORMALIZED TO SUM=1.
                if(MAXC>100) MAXC=100
//              FOR J>100 THERE IS NO COMPOUND NUCLEUS IN PRACTICE.
//              ALEV  IS PARTIAL XSECTION FOR SPIN INDEX I
//              (FROM  COMPOSE  THRU COMMON)
               for(int I=2; I<=MAXC; I++) ALEV[I]+=ALEV[I-1];

                SLEV=ALEV[MAXC];
               for(int I=1; I<=MAXC; I++) ALEV[I]/=SLEV;

//              RUNNING SUM FOR J RANDOM SELECTION FORMED
                JJL=1;
   for(int L=1; L<=NCASC; L++)
            IZCS[L[=IZ;
            INCS[L[=IN;
//              L  IS NO. OF CASCADE BEING CALCULATED
            ERGCS[L]=ENERGY-ECLOSS*RAN(IRX);
            NJL=ALEV[JJL]*NCASC+0.5;
             if(L>NJL)
              do{ JJL++;
                  if(JJL>=MAXC) break
                  NJL=ALEV[JJL]*NCASC+0.5;
                  }while (L>NJL);

             JCS[L]=JJL;
             MJCS[L]=(JJL-1)*MDIR;
          }


//              PRINT 8881,(I,INCS[I],JCS[I],I=1,NCASC)
//          8881FORMAT((4(3I5,5X)))
                goto 230
             210DO 220 L=1,NCASC
                IZCS(L)=IZ
                INCS(L)=IN
                ERGCS(L)=ENERGY
                MJCS(L)=(MAXC-1)*MDIR
             220JCS(L)=MAXC
//              COMPOUND NUCLEUS PREPARATION FINISHED
             230CALL TRACK (LK,LK,LK,EU,EU,EU,ENERGY,1,SPROB,ECLOSS)
                LAT=1
                MAT=NCASC
                IAC=0
                IZC=0
//
//              LOOP TO FIND HIGHEST EX FOR HIGHEST A
//
             240DO 310 I=LAT,MAT
                if(ERGCS[I]==0.) goto 310
                IA=IZCS[I]+INCS[I]
                if(IA-IAC) 310,260,250
             250IAC=IA
                LA=I
//              LOWEST LIMIT FOR A
//              LOOP TO FIND HIGHEST A FOR PARTICULAR Z
                IZC=IZCS[I]
                INC=INCS[I]
                ERGC=0.
                ERGMIN=9999.
//              ERGMIN IS MINIMUM ENERGY FOR GIVEN A
                goto 270
             260IF (IZC==0) goto 250
             270MA=I
//              HIGHEST LIMIT OF A
                IZCSI=IZCS[I]
                INCSI=INCS[I]
                if(IZCSI!=IZC || INCSI!=INC) goto 310
                if(ERGCS[I]-ERGC) 300,290,280
             280ERGC=ERGCS[I]
//              FIND HIGHEST EX FOR HIGHEST A
                LERG=I
             290MERG=I
//              FOR THE ABOVE LA , MA -  LERG , MERG ARE INDICES OF HIGHEST
//              ENERGY LIMITS
             300IF (ERGMIN>ERGCS[I]) ERGMIN=ERGCS[I]
             310}
                if(IAC==0&&LAT==1&&MAT==NCASC) goto 710
                if(IAC<6) goto 710
//
//              WHEN ALL ENERGIES ARE ZERO  -
//              THE WHOLE JOB IS FINISHED
//
                if(IZC==0) goto 230
                SC=.5*MOD(IZC+INC,2)
//              SC IS SPIN INDEX OF NUCLEUS -  0 OR 1/2
               for(320 I=1,4
                II=5-I
                CALL AMASS (RMASS,IZC,INC,II,MEBIN(II),MAXJ(1,II),MAXJS(1,II),EBIN
               1(1,II),RLEV(1,II),IXPR(II),BE)
             320}
                VQS(1)=0.
                VQS(2)=1.44*IZR(2)/((IZR(2)+INR(2))**.3333)
                VQS[3]=2.88*IZR[3]/((IZR[3]+INR[3])**.3333)
//              GET MASS EXCESSES OF NUCLEI
//
//           RM ASS  ARE MASSES OF RESIDUAL NUCLIDES
//          THE EMIN ARE THE SMALLEST ENERGIES FOR WHICH RELIABLE CROSS SECTIONS
//              ARE CALCULABLE.  THE FORMULAE WERE OBTAINED FROM EXPERIENCE.
//
               for(330 I=1,NODES
                EXRSM[I]=ERGC+BE[I]-EMIN[I]
//              EX OF RSID.NUC. - MAXIMUM   MAXIMUM EXCITATION ENERGY OF
//              PARTICULAR RESIDUAL NUCLEUS
             330}
                EXRSM[4]=ERGC
               for(340 K=1,MODES
                if(EXRSM(K)<=0.) goto 340
//              FOR  EX<0.  MODE NOT AVAILABLE
                if(EXRSM(K)<=EBIN(MEBIN(K),K)) goto 340
                WRITE(10,1020)IZR(K),INR(K),EXRSM(K),EBIN(MEBIN(K),K),MEBIN(K)
                STOP
             340}
             350JC=JCS(LERG)
//            J C  IS CASCADE SPIN INDEX
               for(360 I=1,NODES
                if(EXRSM[I]>0.) goto 390
             360}
//
//          EVA PORATION OF PARTICLES OF PARTICULAR Z,N WITH ENERGY E OR LESS NO
//              LONGER POSSIBLE
//
                if(EXRSM[4]<=0.) goto 370
                IBOUND=1
                goto 390
             370LAT=LA
                MAT=MA
               for(380 I=LAT,MAT
//              ZERO THE ENERGY OF THE CASCADE TO AVOID FURTHER CALCULATIONS
                IZCSI=IZCS[I]
                INCSI=IZCS[I]
                if(IZCSI==IZC&&INCSI==INC) ERGCSI=0.
             380}
                IZC=0
                goto 240
//
//              NEW Z,N OF MAX A NEEDED
//
             390DO 410 K=1,MODES
//              A AND EX IN PREDETERMINED LIMITS. J FIXED. PARTICLES
//              AND GAMMAS CAN BE EMITTED.
                if(EXRSM(K)>0.) goto 400
                IXRSM(K)=0
                goto 410
             400CALL FIND (EXRSM(K),IXRSM(K),MAXJ(1,K),MAXJS(1,K),EBIN(1,K),MEBIN
               1(K),K)
//              THE FOLLOWING IS A PROCEEDURE TO FIND THE INTEGER INDEX
//              ASSOCIATED WITH A GIVEN EXCITATION ENERGY IN THE LEVEL
//              FILE. EBIN IS THE ENERGY OF THE BIN IXRSM.
                EXMIN=EXRSM(K)-EMAX(K)+EMIN(K)-1.
                CALL FIND (EXMIN,IXMIN(K),MAXJ(1,K),MAXJS(1,K),EBIN(1,K),MEBIN(K)
               1,K)
             410}
             420DO 430 J=1,2222
//              ZERO THE PROBABILITIES BEFORE CALCULATION
             430PROB[J]=0.
                FPROB=0.
                SPROB=0.
                IPROB=0
                IPROBP=0
//              IPROBP COUNTS PARTICLE OUTPUT
               for(440 K=1,MODES
                PROBM(K)=0.
//              PROBM IS MAXIMUM PROBABILITY FOR MODE (CHANNEL) K.
                if(EXRSM(K)<=0.) goto 440
                AMR=IZC+INC
//              CHNPRB GIVES PROBABILITY FOR ALL CELLS IN MODE K
                CALL CHNPRB (K,IXRSM(K),JC,SC,SPART(K),RLEV(1,K),PROBM(K),MAXJ(1,K
               1),MAXJS(1,K),EBIN(1,K),AMR,IXMIN(K),MEBIN(K),IXPR(K))
             440}
                if(IPROB==0) goto 480
//            A N EXCITED TRAP
                FPROB=PROB(IPROB)*LFISS
               for(450 I=2,IPROB
                PROB[I]=PROB[I]+PROB(I-1)
             450}
                SPROB=PROB(IPROB)
                if(SPROB!=0.) goto 460
                IPROB=0
                goto 480
             460DO 470 I=1,IPROB
             470PROB[I]=PROB[I]/SPROB
                FPROB=FPROB/SPROB
             480MZ=LERG
                MG=MERG
                LERG=0
                IL=0
               for(500 I=MZ,MG
                QIF=Ifabs(INCS[I]-INC)+Ifabs(IZCS[I]-IZC)+fabs(ERGCS[I]-ERGC)
                if(QIF!=0.) goto 500
                if(JC!=JCS[I]) goto 490
                IL=IL+1
                IMGZ(IL)=I
                goto 500
             490IF (LERG==0) LERG=I
             500}
//              ABOVE LOOP SETS UP FOLLOWING:
//              IL   - NUMBER OF CASCADES HAVING RIGHT Z,N,SPIN,EX TO USE SAME
//              PROBABILITY DISTRIBUTION
//              IMGZ - THE RUNNING INDEX NUMBER OF THESE CASCADES
//              LERG - FIRST CASCADE HAVING DIFFERENT SPIN (SAME Z,N,EX)
                if(IL==0) goto 630
               for(620 ICSC=1,IL
                I=IMGZ(ICSC)
//              SCAN ALL EVENTS LOOKING FOR ONES WITH RIGHT A,Z,EX,J
                IZCSI=IZCS[I]
                INCSI=INCS[I]
                JCSI=JCS[I]
                MJCSI=MJCS[I]
//              FOLLOWING INSERT FORCES DECAY THRU E2 GAMMA YRAST CASCADE
                if(IPROB!=0) goto 540
                ERGNOT=ERGCS[I]
             510IF (JCS[I]>1) goto 520
                ERGCS[I]=0.
                JCS[I]=0.
                JCSI=0
                goto 610
             520AJN=JCS[I]-1.+SC
                MODE=4
                ETRAN=ERGCS[I]*(4.*AJN-2.)/(AJN*(AJN+1.))
                ERGC=ERGCS[I]
                ERGCS[I]=ERGCS[I]-ETRAN
                if(ERGCS[I]<0.) ERGCS[I]=0.
                if(JCS[I]!=JC) goto 530
                INOTGS=INOTGS+1
                ERGNOT=ERGC
             530JCSI=JCS[I]
                JCS[I]=JCS[I]-2
//              NOT DECAYED TO G.S. (''YRAST TRAPS''- NO E3 OR M3 IN PROGRAM)
                goto 560
//             }//end YRAST CASCADE INSERT.  15 MARCH 1979
             540CALL SEARCH (IZCS[I],INCS[I],ERGCS[I],JCS[I],IZC,INC,MODE)
                if(ERGCS[I]>=0.) goto 560
                ERGCS[I]=0.
                ERGCSI=ERGCS[I]
                CALL TRACK (5,JCSI,JCSI,ERGC,ERGCSI,0.,ERGC,2,SPROB,ECLOSS)
                if(IDUMP==1) goto 550
                IBUF2(1)=I
                IBUF2(2)=6
//              MODE=6 IS FISSION TRACEBACK
                IBUF2[3]=JCSI
                IBUF2[4]=0
                IBUF2[5]=MJCSI
                IBUF2(6)=ERGC+1.
                IBUF2(7)=0
                IBUF2(8)=IZC
                IBUF2(9)=INC
                WRITE (2) IBUF2,FPROB
             550IZCS[I]=-IZCS[I]
//              NEGATIVE IZCS[I] IMPLIES FISSION OF IZCS[I]
                goto 620
//              SEARCH DETERMINES FINAL NUCLEUS
//
//              PROCEED TO CATALOGUE DISTRIBUTIONS.
//
             560EP=ERGC-ERGCS[I]+BE(MODE)
                if(IPROB==0) EP=ETRAN
                JCSF=JCS[I]
                ERGCSI=ERGCS[I]
                CALL TRACK (MODE,JCSI,JCSF,ERGC,ERGCSI,EP,ERGC,2,SPROB,ECLOSS)
                AC=IZCSI+INCSI
                AP=IZPART(MODE)+INPART(MODE)
                EP=EP*AC/(AP+AC)
                if(IDUMP==1) goto 590
                IE12=0
                IBUF2(1)=I
                if(MODE!=4) goto 580
                JDDJ=JCSF-JCSI
                JDDJ=Ifabs(JDDJ)
                if(JDDJ!=2) goto 570
                IE12=1
                goto 580
             570GMSUM=GE1+GE2*EP**2
                GMSUM=GE1/GMSUM
//              IF DELTA-J IS 2 DECAY IS E2
//              OTHERWISE, DECAY MODE (FOR TRACEBACK PURPOSES ONLY)
//              IS DECIDED BY SIZE OF DECAY PROBABILITY (E1 OR E2 ONLY)
                X=RAN(IRX)
                if(X>GMSUM) IE12=1
             580IBUF2(2)=MODE+IE12
//              IBUF2(2) NOW GIVES MODE=4 FOR E1 GAMMAS AND MODE=5 FOR E2 GAMMAS
//              *** THIS IS FOR THE DUMP-TRACEBACK ONLY
                IBUF2[3]=JCSI
                IBUF2[4]=JCSF
                IBUF2[5]=MJCSI
                IBUF2(7)=EP*10.+1.
                IBUF2(6)=ERGC+1.
                IBUF2(8)=IZC
                IBUF2(9)=INC
                WRITE (2) IBUF2,FPROB
             590JKQ=EP+1.
                if(JKQ>50) JKQ=50
                NSPC(MODE,JKQ)=NSPC(MODE,JKQ)+1
                if(IDIST==0) goto 610
                if(MODE!=MODES) goto 600
                JCSFF=JCSF-1
                if(MJCS[I]<-JCSFF) MJCS[I]=-JCSFF
                if(MJCS[I]>JCSFF) MJCS[I]=JCSFF
                IG=2
                if(IBOUND==1) IG=4
//              STORE GAMMA ENERGY
                CALL OUTEM (IG,MODE,EP,DUM)
                goto 610
             600CALL MJRAN (JCSI,JCSF,MJCSI,MJCSF,MODE,EP,IZCSI,INCSI)
                MJCS[I]=MJCSF
                ID=2
                CALL MOMENT (I,AC,AP,MODE,EP,NCASC,ID)
//              STORE PARTICLE ENERGY. FORMS OUTPUT DISTRIBUTIONS BY
//              CALL ING OUTEM.
             610IF (ERGCS[I]>0.&&IPROB==0) goto 510
             620}
             630CALL TIMEND (NCASC,ISTOP,NTIME)
                if(ISTOP==0) goto 710
                if(IPROB!=0) goto 640
//              PRINT3470,IZC,INC,IAC,ERGC,JC,IBOUND,INOTGS
                ENOTGS=ENOTGS+ERGNOT*(INOTGS-INOTG1)
                SNOTGS=SNOTGS+JC*(INOTGS-INOTG1)
                INOTG1=INOTGS
//              INOTGS=0
//              REMOVE PREVIOUS 2 C-S IN ORDER TO OBTAIN DETAILS ABOUT TRAPS
             640IF (LERG==0) goto 650
//
//              NO FURTHER SPIN AT THIS A AND Z
//              NEW ENERGY OR NUCLIDE NEEDED
//
                JC=JCS(LERG)
                goto 420
//
//              NEW SPIN NEEDED
//
             650IBOUND=0
                LAT=LA
                MAT=MA
                ERGC=0.
               for(680 I=LAT,MAT
                if(ERGCS[I]==0.) goto 680
                if(INC!=INCS[I] || IZC!=IZCS[I]) goto 680
                if(ERGCS[I]-ERGC) 680,670,660
             660ERGC=ERGCS[I]
                LERG=I
             670MERG=I
             680}
                if(ERGC==0.) goto 700
               for(690 K=1,NODES
             690EXRSM(K)=ERGC+BE(K)-EMIN(K)
                EXRSM[4]=ERGC
                goto 350
//
//              SAME Z,N NEW ENERGY NEEDED
//
             700IZC=0
                goto 240
//
//              NEW NUCLIDE OF SAME A NEEDED
//
             710}
                if(IDIST==0) goto 730
               for(720 I=1,NCASC
                AC=IZCS[I]+INCS[I]
             720IF (IZCS[I]>0) CALL MOMENT (I,AC,0.,0,0.,0,3)
//              LAST PARAMETER IN MOMENT
//              = 1  INITIALIZATION OF TABLES
//              = 2  COLLECTING DATA DURING RUN
//              = 3  TERMINATION AND OUTPUT.
             730WRITE(10,1030)
                CALL PRODCT (IDIST,IZMAX,INMAX,NCASC)
                if(IDIST!=0) CALL OUTEM (3,IDUM,DUM,DUM)
                if(INOTGS==0) goto 740
                ENOTGS=ENOTGS/INOTGS
                PR=(100.*INOTGS)/NCASC
                SNOTGS=SNOTGS/INOTGS
                WRITE(10, 1040) PR,ENOTGS,SNOTGS
             740WRITE(10, 1050)
               for(750 IST=1,4
             750BE(IST)=0.
               for(760 IKQ=1,50
                IKQ1=IKQ-1
                IKZ=IKQ
                if(IKQ==50) IKZ=99
                if(NSPC(1,IKQ)!=0 || NSPC(2,IKQ)!=0 || NSPC(4,IKQ)!=0.OR
               1.NSPC(3,IKQ)!=0)WRITE(10,1060)IKQ1,IKZ,(NSPC(JKQ,IKQ),JKQ=1,4)
                EPART=IKQ-.5
               for(760 IST=1,4
                NSPC(IST,51)=NSPC(IST,51)+NSPC(IST,IKQ)
                BE(IST)=BE(IST)+EPART*NSPC(IST,IKQ)
             760}
               for(770 IST=1,4
             770BE(IST)=BE(IST)/(NSPC(IST,51)+1.E-9)
                WRITE(10, 1070) (NSPC(JKQ,51),JKQ=1,4)
                WRITE(10, 1080) BE
               for(780 IST=1,4
                BE(IST)=NSPC(IST,51)
             780BE(IST)=BE(IST)/NCASC
                WRITE(10,1090) BE
                CALL TRACK (LK,LK,LK,ERGC,ERGC,ERGC,ERGC,3,SPROB,ECLOSS)
                WRITE(10, 1100)
                if(ISTOP==0) STOP
//          IBM CALL JSTIME(JTIME2)
//          TMU CALL GETTR (TIME2)
                TT2=SECNDS(TT1)
                CALL LIB$STAT_TIMER(2,ICPU1)
//              CALL CLOCK(JTIME2)
//              T1=(JTIME2-JTIME1)*1.0
//              WRITE(10, 1110) T1
                if(IDUMP==1) STOP
                CALL LIB$INIT_TIMER
                CALL STATIS
//          IBM CALL JSTIME(JTIME3)
//          TMU CALL GETTR (TIME3)
//              CALL CLOCK(JTIME3)
                TT3=SECNDS(TT1)
                T3=TT3-TT2
//              T2=(JTIME3-JTIME2)*1.0
                CALL LIB$STAT_TIMER(2,ICPU2)
                TCPU1=ICPU1/100.
                TCPU2=ICPU2/100.
                WRITE(10,1120) TT2,TCPU1,T3,TCPU2
                STOP
//
  8           10PRINT 9004,LISFIL
                STOP
  8            9PRINT 9006,INFIL
                STOP
  8            2PRINT 9004,EVFIL
                STOP
  8            3PRINT 9005
                STOP
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
*/
void TLLL (double &ER, int NP, int IZ, int IN, double *TL, int &LLMAX)
{
//              WRITTEN BY A. GAVRON
double TLPN[3][17][31], TLAL[31][31], TLT[37][81], DE[4], DE1[4];
int LC[37],LL[4][31];
//char *NMODE[3]={"NEUT","PROT","ALPH"};

//              TRANSMISSION COEFFICIENTS WRITTEN TO TAPE 4
double PROB, SUM, AJ, XSEC, EE, DEE, DTT, T0;
int K=0, IP=1;
double TLL=log(1.E-28);
int IA;
int LC1, NE, LL1, LMT;

if(IP-3<0) {
           IA=IZ+IN;
           K=1;
           DE[IP]=(EMAX[IP]-EMIN[IP])/29.;
           for(int IEN; IEN<=30; IEN++)  EP[IEN]=EMIN[IP]+(IEN-1.)*DE[IP];

// CHECK!!!           TCCAL(IP,IZ,IA,30,15,IMT,ATT,RIV,AIV,AV,AW,AWV,RCX,EP,ETHRS,
// CHECK!!!                        ECUT,ASO,RSOT,ASOT,IRAD,IMESH,TLT,LC);

           for(int IEN=1; IEN<=30; IEN++ ) {
                LC1=LC[IEN];
                SUM=0;
                for(int IJ=1; IJ<=LC1; IJ++) {
                        AJ=IJ-1;
                        SUM+=(2.*AJ+1.)*TLT[IEN][IJ];
                        }

                XSEC=679.58/(EP[IEN]*IA/(IA+1.))*SUM;
                fprintf(f04,"\nIEN=%2d %6.2f  %8.2f", IEN,  EP[IEN],XSEC);
                for(int IJ=1; IJ<=LC1; IJ++) fprintf(f04," %10.2e", TLT[IEN][IJ]);
                LL[IP][IEN]=min(15,LC1-1);
                for(int L=1; L<=16; L++) {
                        TLPN[IP][L][IEN]=log(1.E-28);
                        if( L<=LC1 && TLT[IEN][L]>0)
                                TLPN[IP][L][IEN]=log(TLT[IEN][L]+1.E-28);
                        }
                }
              IP++;
              }
//---------------------------------------------
else if(IP-3==0) {
               DE[3]=(EMAX[3]-EMIN[3])/29.;
               for(int IEN=1; IEN<=30; IEN++) EP[IEN]=EMIN[3]+(IEN-1.)*DE[3];

// CHECK!!!              TCCAL (IP,IZ,IA,30,29,IMT,ATT,RIV,AIV,AV,AW,AWV,RCX,EP,ETHRS,
// CHECK!!!                        ECUT,ASO,RSOT,ASOT,IRAD,IMESH,TLT,LC);

               for(int IEN=1; IEN<=30; IEN++) {
                        LC1=LC[IEN];
                        SUM=0;
                        for(int IJ=1; IJ<=LC1; IJ++) {
                                AJ=IJ-1;
                                SUM+=(2.*AJ+1.)*TLT[IEN][IJ];
                                }
                        XSEC=679.58/(EP[IEN]*4.*IA/(IA+4.))*SUM;
                fprintf(f04,"\nIEN=%2d %6.2f  %8.2f", IEN,  EP[IEN],XSEC);
                for(int IJ=1; IJ<=LC1; IJ++) fprintf(f04," %10.2e", TLT[IEN][IJ]);

                        LL[3][IEN]=min(29,LC[IEN]-1);
                        for(int L=1; L<=30; L++) {
                                TLAL[L][IEN]=log(1.E-28);
                                if(L<=LC1 && TLT[IEN][L]>0) TLAL[L][IEN]=log(TLT[IEN][L]+1.E-28);
                                }
                        }
                for(int IP=1; IP<=3; IP++) DE1[IP]=1./DE[IP];
                VQ0[1]=0;
                VQ0[2]=1.44*(IZ-1.)/pow_dr13(IA-1.);
                VQ0[3]=2.88*(IZ-2.)/pow_dr13(IA-4.);
                IP=4;
                }
//---------------------------------------------
else            {
                EE=ER+VQ0[NP]-VQS[NP]+clshft[NP];
//              THIS PROCEEDURE ADJUSTS THE T-LS TO THE CHANGING COULOMB
//              BARRIER DURING DEEXCITATION.
                if(EE>=EMAX[NP]) EE=EMAX[NP]-.01;
                if(EE<=EMIN[NP]) EE=EMIN[NP]+.01;
                NE=(EE-EMIN[NP])*DE1[NP];
                DEE=EE-EMIN[NP]-DE[NP]*double(NE);
                NE++;
                LLMAX=LL[NP][NE];
                LL1=LLMAX+1;
                if(NP!=3) {
//              LOOPS BELOW REWRITTEN FOR VECTORIZATION
                    for(int LP1=1; LP1<=LL1; LP1++) {
                        TL[LP1]=0;
                        DTT=0;;
                        T0=TLPN[NP][LP1][NE];
                        DTT=(TLPN[NP][LP1][NE+1]-T0)*DE1[NP];
                        TLL=T0+DTT*DEE;
                        TL[LP1]=exp(TLL);
                        }
                }
                else {
                  TL[1]=0;
                  LMT=LLMAX+1;
                  for(int LP1=1; LP1<=LMT; LP1++)  {
                        TL[LP1]=0;
                        DTT=0;
                        T0=TLAL[LP1][NE];
                        DTT=(TLAL[LP1][NE+1]-T0)*DE1[NP];
                        TLL=T0+DTT*DEE;
                        TL[LP1]=exp(TLL);
                        }
                   K++;
                   }
                }
}
/*
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
                double SEARCH (IZ,IN,ERG,J,IZC,INC,MODE)
//
//              ORIGINAL JULIAN double
//
//              ***********
//
//              double  DETERMINES FINAL RESIDUAL NUCLEUS
//
                REAL*8 PROB
                R=RAN(IRX)
//              THE SEARCH ALGORITHM IS TO KEEP HALVING THE VECTOR AND
//              DECREASING LIMITS TILL INDEX IS FOUND. THE INDEX OF THE
//              DETERMINES THE DECAY MODE
                if(R<=PROB(IPROB)) goto 10
                I=IPROB
                goto 70
              10IF (R>PROB(1)) goto 20
                I=1
                goto 70
              20IS=0
                JF=IPROB
              30JF=JF/2+1
                if(JF==2) JF=1
                IS=IS+JF
                if(IS>IPROB) IS=IPROB
              40IF (R>PROB(IS)) goto 30
              50JF=JF/2+1
                if(JF==2) JF=1
                IS=IS-JF
                if(IS<1) IS=1
                if(JF==1) goto 60
                goto 40
              60IF (R<=PROB(IS)) goto 50
                I=IS+1
              70MODE=MJ[I]/19696+1
                if(MODE<5) goto 80
                ERG=-1.
                return
//              IBIN IS PLACE IN 19696 LONG VECTOR OF ONE SPECIFIC MODE.
              80IZ=IZC-IZPART(MODE)
                IN=INC-INPART(MODE)
                MJI=MJ[I]
                IBIN=MOD(MJI,19696)
                if(IBIN==0) IBIN=19696
                MEB=MEBIN(MODE)
               for(90 I=1,MEB
                if(MAXJS(I,MODE)>=IBIN) goto 100
              90}
                WRITE(10, 110) MODE,MEB,MAXJS(MEB,MODE),IBIN
                STOP
             100J=MAXJ(I,MODE)-MAXJS(I,MODE)+IBIN
                ERG=EBIN(I,MODE)
                return
//
             110FORMAT (' SEARCH ERROR.'/' MODE,MEBIN,MAXJS,IBIN'/1X,4I5)
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
                double AMASS (A,IZC,INC,MODE,MEBIN,MAXJ,MAXJS,EBIN,RLEV,IXPR
               1,BE)
//
//              MODIFIED JULIAN double
//
//              **************
//              DETERMINE MASSES.
                DIMENSION MAXJ(286), MAXJS(286), EBIN(286), RLEV(19696), BE(4), A(
               14)
                MAXXL=MAXC+9
                if(MAXXL>110) MAXXL=110
                IZR(MODE)=IZC-IZPART(MODE)
                INR(MODE)=INC-INPART(MODE)
                CALL MASSES (IZR(MODE),INR(MODE),A(MODE))
                BE(MODE)=A[4]-A(MODE)-PMASS(MODE)
                if(MODE!=4) goto 20
                WRITE (4,30)
                IA=IZR(MODE)+INR(MODE)
                IZ=IZR(MODE)
               for(10 I=1,MAXXL
                L=I-1
                CALL BARFIT (IZ,IA,L,FBARR[I],EROT[I],SELMAX)
                FBARR[I]=FBARR[I]*BARFAC
                EROT[I]=EROT[I]*FYRST
              10}
              20CALL RANGE (MODE,IMIN,IMAX)
//        ccccc modified by ikezoe 1988 7/23
                aratio10=1.
                if(izr(mode)==izshel.and.inr(mode)==inshel)then
                aratio10=pshel
               }//endif
                CALL LEVDEN (IZR(MODE),INR(MODE),IMIN,IMAX,MAXJ,MAXJS,EBIN,MEBIN
               1,RLEV,IXPR,BE(MODE),aratio10)
                return
//
              30FORMAT (1H )
               }//end */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PRODCT (int IDIST, int IZ, int IN, int NCASC)
{
//
//              ORIGINAL JULIAN double
//
//              ***********
//
//              SORTING ROUTINE OF RESIDUAL NUCLEI
//
//------COMMON /OUT/ -------------------
/*
      COMMON /OUT/ IZCS(9996), INCS(9996), JCS(9996), ERGCS(9996), MJCS(9996), RLEV(19696,4)
      COMMON /XQSIG/ SIGMA
      COMMON /XQANG/ SUM(300,10), MDIR, CT
*/
int NRES[991], NNUC[991], IZNUC[991], INNUC[991];
double PRES[991];
double SIG;
int ID, JK, IK, I, J;

int *KNUC=MJCS;

               for(J=1; J<=990; J++ ) NRES[J]=0;

                int IRES=0;
                int IFISS=0;
                int NC=NCASC;
               for(I=1; I<=NC; I++ ) {
                      if(IZCS[I]<0)  IFISS++;
                      else  {
                        int IZC=IZ-IZCS[I];
                        int INC=IN-INCS[I];
                        int IA=IZC+INC;
                        int IBIN=IA*(IA+1)/2+IZC+1;
                        if(IBIN>990) {
                              if(IRES!=1) {
                                fprintf(f10,"\nTHE FOLLOWING PRODUCTS WERE NOT BINNED   Z    N");
                                IRES=1;
                                }
                              else fprintf(f10,"\n %3d %5d", IZCS[I],INCS[I]);
                              }
                        else NRES[IBIN]++;
                      }
                   }

                fputs(s100m,f10);
                fprintf(f10,"\nYields of residual nuclei %s    Z      N     X     A   "
                            "Events      Percent      X-section(mb)",s35x);
                fputs(s100m,f10);
                double SRES=NCASC;
                IAT=IZ+IN;
                int K=1;
                double SIGS=0;
                int NRESS=0;
                double PRESS=0;

           for(I=1; I<=44; I++ ) { //   IF DIMENSION CHANGED FROM 990, CHANGE THIS      for(RANGE.
                int IBIN1=I*(I-1)/2+1;
                int IBIN2=I*(I+1)/2;
                int IZR=IZ;
                IAT=IZ+IN-I+1;
                for(J=IBIN1; J<=IBIN2; J++) {
                        if(NRES[J]!=0) {
                                PRES[J]=NRES[J]/SRES*100.;
                                int INR=IAT-IZR;
                                SIG=SIGMA*PRES[J]/100.;
                                fprintf(f10,"\n%s%3d   %3d   %3d %.2s  %6d  %7.3f  %9.2f",
                                         s35x,IZR,INR,IAT,LMNT(IZR),NRES[J],PRES[J],SIG);

                                if(PRES[J]>=3.) {
                                        IZNUC[K]=IZR;
                                        INNUC[K]=INR;
                                        NNUC [K]=NRES[J];
                                        K++;
                                        }
                                SIGS+=SIG;
                                PRESS+=PRES[J];
                                NRESS+=NRES[J];
                                }
                        IZR--;
                        }
                }

                double PR=IFISS/SRES*100.;
                SIG=SIGMA*PR/100.;
                PRESS+=PR;
                NRESS+=IFISS;

                if(IFISS>0)
                    fprintf(f10,"\n%sTotal fission      %6d   %7.3f  %9.2f",s35x,IFISS,PR,SIG);

                SIGS+=SIG;
                fprintf(f10,"\n%s---------------------------------------",s35x);
                fprintf(f10,"\n%s       TOTAL:   %6d    %7.3f   %9.2f", s35x,NRESS,PRESS,SIGS);
                fputs(s100m,f10);

                if(IDIST==0) return;
                ID=1;
                int IDUM=0;                       // check
                OUTRES (ID,IDUM,IDUM);
                K--;
               for(I=1; I<=NC; I++) {
                     if(IZCS[I]<0) continue;

                     JK=16;
                     for(J=1; J<=K; J++)
                           if(IZCS[I]==IZNUC[J] && INCS[I]==INNUC[J])  {KNUC[I]=JK=J;}
                     IK=I;
                     ID=2;
                     OUTRES (ID,IK,JK);
                     }


                K++;
                fprintf(f10,"\n******  ANGULAR DISTRIBUTION RESULTS ******"
                                "************************************");
                K=min(16,K);
                if(MDIR==0) fprintf(f10,"\n\n***** Spin alignment perpendicular to recoil axis\n"
                                "- standard compound nucleus angular distribution");

                if(MDIR==1) fprintf(f10,"\n\n***** Spin alignment perpendicular to reaction plane\n"
                        "- angular distribution is around Z-axis perpendicular to reaction plane");

               for(I=1; I<=K; I++) {
                       fputs(s100m,f10);
                       if(I!=K)  fprintf(f10,"\n%sEnergy and angular distribution "
                              "of residual nucleus - Z=%d  N=%d",s35x,IZNUC[I],INNUC[I]);
                       else    fprintf(f10,"\n%sEnergy and angluar distribution "
                              "of all residual nuclei");

                        fputs(s100m,f10);
                        ID=3;
                        JK=J;
                        IK=I;
                        if(IK==K) IK=16;
                        OUTRES (ID,JK,IK);
                        }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void COMPND (int IZ, int IN, double ENERGY, int &MAXC)
{
//
//              ORIGINAL JULIAN double
//
//              ***********
//              FORM J DISTRIBUTION OF INITIAL COMPOUND NUCLEI
extern void GCLDP (double AMR, double AZR, int IADEF, double FALIT, double &ALIT,
            double &SIGSQ, double FACT1, double& DELTA, double &ACONST, double SR, double ARATIO);

double t;
double Z=IZ;
double A=IZ+IN;

double SR=double((IZ+IN)%2)*0.5;
double FL=0;   // check FL
double AL=0;   // check AL
double SIG=0;  // check SIG
double F1=0;   // check F1
double D =0;   // check D
double ACONST=0;   // check ACONST

GCLDP(A,Z,0,FL,AL,SIG,F1,D,ACONST,SR,1.);

double UEX=ENERGY-D;
 SIG=1./(2.*SIG*sqrt(UEX/AL));

MAXC=min(MAXC,100);

for(int J=1; J<=MAXC; J++) {
              double AJ=J-1.+SR;
              t = AJ+0.5;
              ALEV[J]=(2.*AJ+1.)*exp(-t*t*SIG);
              }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void COMPOS(int &IZ,int &IN,double &EEXCN,double &VZC,int &MAXC,
                double &SIGCN,double &ECLOSS,int INPUTp)
{
//           AMP     IS A OF PROJECTILE
//           AZP     IS Z OF PROJECTILE
//           SP is SPIN OF PROJECTILE
//           ST is SPIN OF TARGET
//           QCN     IS Q OF COMPOUND NUCLEUS FORMATION
//           IPOT    IS OPTICAL MODEL POTENTIAL CODE

double TLCN[151], CSCH[9];
int INP, INT,IZCM,INCM, NTEMP, I;
//              **************************

if(INPUTp<=1) {
                        //     READ (9,410) IZP,IAP,IZT,IAT,SP,ST,QCN
     if(QCN==0) {
             INP=IAP-IZP;
             INT=IAT-IZT;
             IZCM=IZP+IZT;
             INCM=INP+INT;
             NTEMP=NOSHL;
             NOSHL=0;
             double DFCP=MASSES (IZP,INP);
             double DFCT=MASSES (IZT,INT);
             NOSHL=NTEMP;
                //  MODIFICATION 1/31/83  G.S. OF PROJECTILE AND TARGET ARE ALWAYS SHELL CORRECTED.
             double DFCCM=MASSES (IZCM,INCM);
             QCN=DFCP+DFCT-DFCCM;
             }

 int  IDENT=0;
   if(IZP==IZT && IAP==IAT ) IDENT=1;
 int  ISPIN=0;
   if(SP!=0. || ST!=0.) ISPIN=1;
// !!! CHECK OLEG   if(IDENT==1&&ISPIN==1)WRITE(10,400)
//              **************************

   double  ELB1=ELAB-0.5*ELOSS;
   int     IXPER=EXPSIG;

   if(EXPSIG==0.) EXPSIG = max(BASS(IAP,IZP,IAT,IZT,ELB1),1);
//              AGRAZ IS DIFFUSENESS OF TRANSMISSION COEFFICIENTS.
//              IF AGRAZ!=0  TLOM IS BYPASSSED
//              EXPSIG IS CORRECT EXPERIMENTAL XSECTION IF EXISTS. IT IS USED
//              TO CUT THE TL DISTRIBUTION TO PROVIDE THE CORRECT XSECTION.

   double  AZP=IZP;
   double  AMP=IAP;
   double  AZT=IZT;
   double  AMT=IAT;
           IZ=IZP+IZT;
   int     INP=IAP-IZP;
   int     INT=IAT-IZT;
           IN=INP+INT;
   int     IA=IZ+IN;
   int     AMC=IA;
   double  ECM=ELAB*AMT/(AMP+AMT);
   double  AMU=AMP*AMT/(AMP+AMT);
   double  ALAMDA=4.651/sqrt(AMU*ECM);
           ECLOSS=ELOSS*AMT/(AMP+AMT);
   double  BRR,SEGS,SELMAX;
   BARFIT (IZ,IA,0,BRR,SEGS,SELMAX);

   double t=ALAMDA*(SELMAX+1.);
   double SIGMAT=31.4159*t*t;

//              SELMAX IS ANG MOM AT WHICH BARRIER GOES TO ZERO
//              SIGMAT IS MAXIMUM CROSS SECTION FOR WHICH STATISTICAL MODEL IS
//              VALID. HIGHER PARTIAL WAVES HAVE ZERO FISSION BARRIER
   double SIGTMP=EXPSIG;
   int ILIMIT=0;
   if( !(EXPSIG<=SIGMAT || SELMAX<10.) ) {
                EXPSIG=SIGMAT;
                ILIMIT=1;
                }
                                        //  ILIMIT=1 IMPLIES CROSS SECTION  W A S  TRUNCATED
                                        //  ILIMIT=0 IMPLIES CROSS SECTION NOT MODIFIED
   double EC1=ECM-ECLOSS*.5;
   EEXCN=QCN+ECM;
   double ROT=34.548*pow(AMC,-1.6667);
   int LYRAST=sqrt(EEXCN/ROT)+.5;
                                        //  LYRAST IS SPIN OF YRAST (RIGID BODY)
                                        // LEVEL AT C.N. EXCITATION ENERGY.
   AGRAZ=max(0.3,AGRAZ);
   int LL=100;
   int L1=100;
   double EX;

   for(I=1; I<=100; I++) {
             EX=I-40.;
             EX/=AGRAZ;
             if(EX>20.)         TLCN[I]=0;
             else   if(EX<-20.) TLCN[I]=1;
                    else        TLCN[I]=1./(1.+exp(EX));
             }

   double SIGOM=0;
   double AM=AMP*AMT/(AMP+AMT);
   double AK=0.218898*sqrt(AM*EC1);
   double PLSQ=10.*PI/AK/AK;
   double FAC,AI1,SIGP;

   for(I=1; I<=L1; I++) {
                FAC=1;
                AI1=I-1;
                if(IDENT!=0) FAC=(I%2==1  ? 2 : 0 );
                SIGP=(2.*AI1+1.)*TLCN[I]*FAC;
                SIGOM+=SIGP;
                }

   SIGOM*=PLSQ;
//
//           NORMALIZATION OF OPTICAL MODEL AND EXPERIMENTAL CROSS SECTIONS
//
//            T = 1/(1+exp((L-L(G))/V)
//            M*L+B = Y WHERE Y = LN(1/T-1), M = 1/V, AND B = -M * L(G)
//            SOLVE FOR L(G) AND V FROM CALCULATED TRANSMISSION COEFFICIENTS
//            S = C * SUM((2 * L + 1) * T  WHERE L(G) IN T IS NOW L(C)
//           (L(G)+1)**2/(L(C)+1)**2 = SIGMA(CALCULATED)/SIGMA(EXPERIMENTAL)
//            IS USED TO SOLVE FOR L(C)
//
  double  EXPS=EXPSIG/PLSQ;
  double  TLLS=0;
  double  TLS=0;
  int     IJUMP=0;

  int L0,AI;
  int  L2=0;
  double TH;

    for(I=2; I<=L1; I++) {
         if(TLCN[I]>0.9*TLCN[1]) continue;
         if(TLCN[I]<0.004) break;
         L2=I;
         if(IJUMP!=1) {IJUMP=1; L0=I-1;}
         AI=I-1.;
         TH=log(1./TLCN[I]-1.);
         TLS+=TH;
         TLLS+=TH*AI;
         }


    int LS=(L2-1)*L2/2;
    int L2S=LS*(2*L2-1)/3;


    if(L0>1)  {
        int L0S=(L0-1)*L0/2;
        LS-=L0S;
        L2S-=L0S*(2*L0-1)/3;
        }

    double SL=LS;
    double SL2=L2S;
    double AL=L2-L0;
    double DEN=SL*SL-AL*SL2;
    double ELAM=DEN/(TLS*SL-AL*TLLS);
    double ELG=ELAM*(TLS*SL2-SL*TLLS)/DEN;
//    double ELGOM=ELG;
    int ITER=0;

    double SF,EXPON, FUNC, RSIG;

L140: ITER++;
      if(ITER!=101) {//==============
            SF=0.;
            for(I=1; I<=150; I++) { //-------------------
                AI=I-1;
                EXPON=(AI-ELG)/ELAM;

                if(EXPON > 40) break;

                if(EXPON<=-40) DEN=1;
                else           DEN=1.+exp(EXPON);

                FUNC=(2.*AI+1.)/DEN;
                if(FUNC<0.001) break;
                SF+=FUNC;
                }//-----------------------------------

                RSIG=EXPS/SF;
                if(fabs(RSIG-1.)>=0.0001) {
                        ELG=(ELG+1.)*sqrt(RSIG)-1.;
                        goto L140;
                        }
            }        //==============

         L1=0;
         for(I=1; I<=150; I++) {
                AI=I-1;
                EXPON=(AI-ELG)/ELAM;
                EXPON=max(-60,EXPON);
                EXPON=min( 60,EXPON);  //60 DEPENDS ON COMPUTER IN USE.

                TLCN[I]=1./(1.+exp(EXPON));
                if(TLCN[I]<0.001) break;
                L1=I;
                }

            SIGCN=0;
            L1=min(L1,100);


double SCHX=SP+ST;
double SCHN=fabs(SP-ST);
double SIN=(2.0*SP+1.0)*(2.0*ST+1.0);
int    ISC=SP+ST+0.01;
double ASC=ISC;
double SC=SP+ST-ASC;
double AC1, AJC, SUMTL, ALMAX, ALMIN, DL;
int  LMAX,LMIN, L;
double S1,S2,DS;

    //-------------------------
         for(I=1; I<=L1; I++) {
                AC1=I-1;
                AJC=AC1+SC;
                SUMTL=0;
                ALMAX=AJC+SCHX;
                LMAX=ALMAX+0.01;
                if(SCHX<AJC) {
                        ALMIN=AJC-SCHX;
                        LMIN=ALMIN+0.01;
                        }
                else    {
                        if(SCHN>AJC) {
                            ALMIN=SCHN-AJC;
                            LMIN=ALMIN+0.01;
                            }
                        else LMIN=0;
                        }

                L=LMIN;
                do {
                        if(L>LL || L>L1) break;
                        DL=fabs(AJC-L);
                        SL=AJC+L;
                        S1=max(SCHN,DL);
                        S2=min(SCHX,SL);
                        DS=S2-S1+1.;
                        if(DS>=1.0) SUMTL+=DS*TLCN[L+1];
                        L++;
                   } while( L<=LMAX);

                AJCN[I]=AJC;
                SIGJC[I]=(2.*AJC+1.)*PLSQ*SUMTL/SIN;
                if(IDENT!=0) SIGJC[I]=(I%2==1 ?  2*SIGJC[I] : 0 );
                SIGCN+=SIGJC[I];       // PARTIAL XSECTION (IS  ALEV  IN MAIN PROGRAM)
                }
    //-------------------------

double      EREC=ELAB*AMP/AMC;
            VZC=sqrt(2.0*EREC/(931.5*AMC));
double      VZCC=VZC*30.;
double      VZB=sqrt(2.0*ELAB/(931.5*AMP));
double      VZBC=VZB*30.;
            MAXC=L1;
            
char *s450=
"-------------------------------------------------------------------------\n"
"STARTING CONDITIONS\n\n"
"                     Z      N      A     SPIN\n";

fputs(s450,f10);
fprintf(f10,"      Projectile   %3d    %3d    %3d   %5.1f\n",IZP,INP,IAP,SP);
fprintf(f10,"      Target       %3d    %3d    %3d   %5.1f\n",IZT,INT,IAT,ST);
fprintf(f10,"      Compound     %3d    %3d    %3d\n",        IZ,IN,IA);
fprintf(f10,"\n%sBombarding energy                        %10.3f  MeV\n", s35x,ELAB);
fprintf(f10,  "%sCenter of mass energy                    %10.3f  MeV\n", s35x,ECM);
fprintf(f10,  "%sCompound nucleus excitation energy       %10.3f  MeV\n", s35x,EEXCN);

if(ECLOSS>0)
  fprintf(f10,"%sExcitation energy loss thru target       %10.3f  MeV\n", s35x,ECLOSS);
fprintf(f10,  "%sCompound nucleus recoil energy           %10.3f  MeV\n", s35x,EREC);
fprintf(f10,  "%sCompound nucleus recoil velocity         %10.4e  cm/nanosec\n", s35x,VZCC);

fprintf(f10,  "%sCompound nucleus velocity/c              %10.4e\n", s35x,VZC);
fprintf(f10,  "%sBeam velocity                            %10.4e  cm/nanosec\n", s35x,VZBC);
fprintf(f10,  "%sBeam velocity/c                          %10.4e\n", s35x,VZB);

if(EXPSIG!=0.) {

  if(LMINN>0) fprintf(f10,  "%sPartial waves below L = %i exluded\n", s35x,LMINN);
  EXPSIG=SIGTMP;

  if(IXPER>0) fprintf(f10,  "%sExperimental fusion cross section        %10.3f  mb\n", s35x,EXPSIG);
  if(IXPER==0)fprintf(f10,  "%sBass 1977    fusion cross section        %10.3f  mb\n", s35x,EXPSIG);

 if(ILIMIT==1){
              fprintf(f10,  "%sCross section truncated for calculation to %8.2f\n", s35x, SIGMAT);
              fprintf(f10,  "%shigher partial waves have zero fission barrier\n", s35x);
        }
fprintf(f10,  "%sFusion L-grazing                         %10.3f\n", s35x, ELG);
fprintf(f10,  "%sFusion L-diffuseness                     %10.3f\n", s35x, ELAM);
fprintf(f10,  "%sYrast spin at maximum excitation energy  %10d\n", s35x,LYRAST);
fprintf(f10,  "%sCompound nucleus formation cross section %10.3f  mb\n", s35x,SIGCN);
  }

  SC=(IA%2)*0.5;
  if(LMINN==0) return;

  for(LL=1; LL<=LMINN; LL++ ) {
                             //A NG. MOM. LMINN LESS 1 IS INDEX LMINN
       SIGCN-=SIGJC[LL];
       SIGJC[LL]=0;
       }

fprintf(f10,  "%sComp.nuc.cross section in window         %10.3f  mb\n", s35x,SIGCN);
  }/*
 ///------------------------------------------------------------------ part 2
             350READ (9,670) IZ,IA,EEXCN,EREC,AJNUC,LMINN
//              AJNUC IS REAL TRUE SPIN. (NOT THE INDEX)
                IN=IA-IZ
                AMC=IA
                VZC=sqrt(2.0*EREC/(931.5*AMC))
                VZCC=VZC*30.
                SC=MOD(IA,2)*0.5
                AMAXC=AJNUC-SC+1.1
                if(AMAXC<1)AMAXC=1
                MAXC=AMAXC
                if(INPUTp==3) READ (9,680) MAXC,(SIGJC[I],I=1,MAXC)
                if(INPUTp==4) CALL COMPND (IZ,IN,EEXCN,MAXC)
                WRITE(10, 690)
                WRITE(10, 470)
                WRITE(10, 500) IZ,IN,IA
                WRITE(10, 530) EEXCN
                WRITE(10, 550) EREC
                WRITE(10, 560) VZCC
                WRITE(10, 570) VZC
                if(INPUTp>2) goto 360
                WRITE(10, 700) AJNUC
                SIGCN=100.
                return;
             360 SIGCN=0.
                LM1=LMINN+1
               for(370 I=1,MAXC
                if(INPUTp==5) SIGJC[I]=0.
                AC1=I-1
                AJCN[I]=AC1+SC
                if(INPUTp==5&&I>=LM1) SIGJC[I]=2.*AJCN[I]+1.
                SIGCN=SIGCN+SIGJC[I]
             370}
                WRITE(10, 650) SIGCN
             380 return
  */             }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*
                double CHNPRB (IPOT,IXR,ISC,SC,SE,RLEV,PROBM,MAXJ,MAXJS,EBIN
               1,AMR,IXMIN,MEBIN,IXPR)
//
//              MODIFIED JULIAN double
//
//              ***********
                DIMENSION MAXJ(286), MAXJS(286), EBIN(286), RLEV(19696), MAXJT(286
               1), MAXJST(286), EBINT(286), WUE(3), WUM(3), TL(30), TLSUM(30),
               2TLS0(25)
                EQUIVALENCE (TLS0(2),TLSUM(1))
                COMMON /RLVV/ RLEVT(19696)
              
                REAL*8 PROB
//
//           WE ISSKOPF UNITS FOR GAMMA DEEXCITATION
//
                DATA WUE, WUM /6.8E-08,4.9E-14,2.3E-20,2.1E-08,1.5E-14,6.8E-21/,
               1IZPREV, INPREV /0,0/, TLS0 /25*0./
                DATA PI2, TOTHRD /6.28318,.666667/
                SR=fabs(SC-SE)
//              ISC = SPIN OF LEVEL
//              SE = PARTICLE SPIN
//              SC  = BASIC SPIN UNIT  ( 0 OR 1/2 )
//              SR  = MINIMUM SPIN WE CAN DECAY TO
                IXRSM=IXR
//              STARTING FROM HIGHEST LEVEL IXRSM , WE GO LEFT (IGO=-1)
//              OR RIGHT TO LOOK FOR MAXIMUM PROBABILITY.
                ICUT=0
                if(IXMIN<IXPR) ICUT=1
                JSC=ISC
                JMAX=ISC
                MODUL=(IPOT-1)*19696
                AJC=ISC-1+SC
                SEP01=SE+.01
                SEM01=SE-.01
                SRM1=SR-1.
              10IF (IPOT==4) goto 80
                IPROBO=IPROB
                ECM=ERGC+BE(IPOT)-EBIN(IXRSM)
                IZZ=IZR(IPOT)
                INN=INR(IPOT)
                CALL TLLL (ECM,IPOT,IZZ,INN,TL,LLMAX)
                if(TL(1)<1.E-10) goto 70
                TLSUM(1)=TL(1)
                LLM=LLMAX+1
               for(20 I=2,LLM
              20TLSUM[I]=TLSUM(I-1)+TL[I]
                IGO=-1
//              JSR = J OF RESIDUAL NUCLEUS.
                JSR=JSC
                if(JSR>MAXJ(IXRSM)) JSR=MAXJ(IXRSM)
                MDELJ=MAXJS(IXRSM)-MAXJ(IXRSM)
              30AJR=JSR+SRM1
                SUMTL=0.
                DJ=fabs(AJR-AJC)
                SJ=AJR+AJC
                ALMAX=SJ+SEP01
//              ALMAX = SUM OF SPINS
//              LMAX=ALMAX+0.01
                LMAX=ALMAX
                if(LMAX>LLMAX) LMAX=LLMAX
//              MINIMUM L CONTRIBUTING TO TRANSITIONS
                ALMIN=DJ-SEM01
                LMIN=ALMIN
                if(LMIN>LMAX) goto 50
//          +++ +++
//          +++ +++ FOLLWING IS SUBSTITUTE CODE FOR CHANNEL SPIN
//          +++ +++ OF  0  AND  1/2 .
                LM1=LMIN
                LM2=LMAX+1
                SUMTL=TLSUM(LM2)-TLSUM(LM1)
//              FOR LM1=0 TLSUM IS 0 THROUGH EQUIVALENCE TO TLS0
                SUMTL=SUMTL+2.*SE*(SUMTL-TL(LM2)-TL(LM1+1))
//
//              SUMTL IS BEING MULTIPLIED BY  1 FOR SE=0  ALL L-S
//              2 FOR SE=1/2 L ABOVE LMIN BELOW LMAX
//              1 FOR SE=1/2 L=LMIN OR LMAX
//          +++ +++}//end OF SUBSTITUTE CODE.
                MJS=MDELJ+JSR
//              ORDINAL NO. OF THIS ENERGY LEVEL
//          PAR PRO=RLEV(MJS)*SUMTL*.5
//              IGNORE PARITY FACTOR.DOES NOT EFFECT RESULTS
                PRO=RLEV(MJS)*SUMTL
                if(PRO<=DISCPR*PROBM) goto 50
//              SEE REMARK ABOUT DISCPR IN MAIN PROG.
                IPROB=IPROB+1
//              IPROB IS TOTAL NO. OF PROBABILITIES CALCULATED
//              if(IPROB>2222)goto 60
                PROB(IPROB)=PRO
                MJ(IPROB)=MJS+MODUL
                if(PRO<PROBM) goto 40
                PROBM=PRO
                JMAX=JSR
//              SEARCH IN ALL EX AND J DIRECTIONS FOR PROB>=.001PROBMAX
              40IF (JSR==1) goto 60
                JSR=JSR+IGO
                if(JSR>MAXJ(IXRSM)) goto 70
                goto 30
              50IF (IGO==1) goto 70
                if(JSR==JSC) goto 40
              60IGO=1
                JSR=JSC+1
                if(JSR<=MAXJ(IXRSM)) goto 30
              70IF (IXRSM==IXMIN) return
                IXRSM=IXRSM-1
                if(IPROB==IPROBO&&ICUT==0&&TL(1)>=1.E-10) return
                JSC=JMAX
                goto 10
              80AMRT=AMR**TOTHRD
                AMRTT=AMRT**2
                PAMRT=PI2*AMRT
                GE1=FGE1*PAMRT*WUE(1)
                GM1=FGM1*WUM(1)*PI2
                GE2=FGE2*PI2*AMRTT*WUE(2)
                GM2=FGM2*PAMRT*WUM(2)
                ECM=ERGC-EBIN(IXRSM)
                ECUBE=ECM**3
                EFIVE=ECUBE*ECM*ECM
                GGE1=GE1*ECUBE
                GGM1=GM1*ECUBE
                GGE2=GE2*EFIVE
                GGM2=GM2*EFIVE
//          PAR PROBA=(GGE1+GGM1)*.5
//          PAR PROBB=(GGE2+GGM2)*.5
                PROBA=GGE1+GGM1
                PROBB=GGE2+GGM2
//              IGNORE PARITY
                PROBC=PROBA+PROBB
                L=JSC-3
               for(140 I=1,5
                JR=I+L
                if(JR<=0) goto 140
                if(JR>MAXJ(IXRSM)) goto 140
                if(I==1 || I==5) goto 100
                if(I==2 || I==4) goto 110
                if(AJC<.1) goto 140
                if(AJC<.6) goto 120
              90PROB1=PROBC
                goto 130
             100PROB1=PROBB
                goto 130
             110AJR=JR-1+SR
                if(AJC<.1&&AJR<.1) goto 120
                goto 90
             120PROB1=PROBA
             130MJS=MAXJS(IXRSM)+JR-MAXJ(IXRSM)
                PRO=RLEV(MJS)*PROB1
                if(PRO<=DISCPR*PROBM) goto 140
                IPROB=IPROB+1
                if(IPROB>2222) goto 150
                MJ(IPROB)=MJS+MODUL
                PROB(IPROB)=PRO
                if(PRO>PROBM) PROBM=PRO
             140}
                if(IXRSM==IXMIN) goto 160
                IXRSM=IXRSM-1
                JSC=JMAX
                goto 80
             150WRITE(10, 200)
                STOP
             160BARR=FBARR(ISC)
                LFISS=0
                if(ERGC>=BARR+EBIN(1)) goto 161
                IPROB=IPROB+1
                if(IPROB>2222)goto 150
                MJ(IPROB)=IMODL
                PROB(IPROB)=0.
//              BUG FIX 12/31/85
                return
             161IF (IZPREV==IZR[4]&&INPREV==INR[4]) goto 170
                IZPREV=IZR[4]
                INPREV=INR[4]
                IXMAX=ERGC+2.5
                IX=1
                MX=0
                CALL LEVDEN (IZR[4],INR[4],IX,IXMAX,MAXJT,MAXJST,EBINT,MEBINT
               1,RLEVT,MX,BE[4],ARATIO)
             170EXRES=ERGC-BARR
                CALL FIND (EXRES,IXRSM,MAXJT,MAXJST,EBINT,MEBINT,1)
                SUM=0.
               for(180 I=1,IXRSM
                N=IXRSM-I+1
                if(ISC>MAXJT[N]) goto 180
                NLEV=MAXJST[N]-MAXJT[N]+ISC
                R=RLEVT(NLEV)
                SUM=SUM+R
                if(R<.001*SUM) goto 190
             180}
//         P204 PRO=SUM*.5
             190PRO=SUM
//              IGNORE PARITY
                LFISS=1
                IPROB=IPROB+1
                if(IPROB>2222) goto 150
                MJ(IPROB)=IMODL
                PROB(IPROB)=PRO
                return
//
             200FORMAT ('0DIMENSION OF PROB NOT HIGH ENOUGH - ****************
               1')
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

                double MOMENT (I,A,AP,MODE,EP,NCASC,ID)
//
//              ORIGINAL JULIAN double
//
//              ***********
//              MONITORING OF OUTGOING PARTICLES.
                COMMON /MOM/ VY(9996), VZ(9996), VX(9996)
//              EP PARM. COMES IN FROM VZC WHEN ID=1
                DIMENSION ERES(9996), FCT(9996)
                EQUIVALENCE (ERES(1),VY(1)), (FCT(1),VZ(1))
                if(ID-2) 10,30,40

//
//           INITIALIZATION
//
              10N=1
                NN=N+NCASC-1
               for(20 J=N,NN
                VY[J]=EP*MDIR
                VX[J]=0.
              20VZ[J]=EP*(1-MDIR)
                goto 60
//
//           CA LCULATION OF MOMENT
//
              30 RN4=RAN(IRX)
                PHI=6.28318*RN4;
                SOX=2.*AP*EP/(A*A*931.5);
//              BUG FIX  A*(AP+A) REPLACED BY A**2 12/15/82
//              if(SOX<0.)PRINT 943,AP,EP,A
                VT=sqrt(SOX)
                SOX=1.-COSTH*COSTH;
                if(SOX<0.) SOX=0.
                SINTH=sqrt(SOX)
                VZSE=VT*COSTH
                VYSE=VT*SINTH*sin(PHI)
                VXSE=VT*SINTH*cos(PHI)
                VZ[I]=VZ[I]+VZSE
                VY[I]=VY[I]+VYSE
                VX[I]=VX[I]+VXSE
                VZP=VZ[I]-VZSE*(A/AP+1.)
                VYP=VY[I]-VYSE*(A/AP+1.)
                VXP=VX[I]-VXSE*(A/AP+1.)
                VPP2=VZP*VZP+VYP*VYP+VXP*VXP
                EPART=AP*VPP2*469.
                SOX=0.
                if(VPP2>0.0) SOX=VZP/sqrt(VPP2)
//             for(NOT USE QUICK FUNCTIONS HERE
                ANG=Acos(SOX)*180./3.1415927
                CALL OUTEM (2,MODE,EPART,ANG)
                goto 60
//
//           EN D CALCULATION
//
              40 VFTS=VX[I]*VX[I]+VY[I]*VY[I]+VZ[I]*VZ[I];
                ERES[I]=0.5*A*VFTS*931.5;
                if(VFTS==0.)   FCT[I]=0.0
                else           FCT[I]=cos(VZ[I]/sqrt(VFTS))*180./PI;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
*/
/*
COMMON /XQSIG/ SIGMA, NCASC, DSIG(18)
*/

void OUTEM (int ICTL, int MODE, double EEML, double AEML)
{
//
//
//           ICTL IS ENTRY CONTROL
//           MODE IS KIND OF PARTICLE EMITTED OR GAMMA
//           EEML IS LABORATORY ENERGY OF EMITTED PARTICLE OR GAMMA
//           AEML IS LABORATORY ANGLE OF EMITTED PARTICLE
//           NTIS NUMBER OF EVENTS OF EACH MODE
//           N( MODE,ENERGY,ANGLE) IS DISTRIBUTION OF EVENTS
//              61ST ENERGY IS FOR ENERGY OVERFLOW
//              62ND ENERGY IS FOR TOTAL AT EACH ANGLE
//              THERE ARE 18 ANGLES
//           NG(ENERGY) IS FOR EACH ENERGY NUMBER OF GAMMAS
//           NW(MODE,BIN,ANGLE) IS DISTRIBUTION OF EVENTS INTO FOUR ENERGY BINS
//           NGW(BIN) IS DISTRIBUTION OF GAMMAS INTO FOUR ENERGY BINS
//           EW IS ENERGY OF SEPARATION BETWEEN FOUR ENERGY BINS
//           EW1 IS LOWER ENERGY OF EACH OF FOUR ENERGY BINS
//           EW2 IS UPPER ENERGY OF EACH OF THREE ENERGY BINS
//              THE FOURTH IS INFINITY
//
//           DELE IS THE ENERGY DIFFERENCE FOR THE 30 ENERGY BINS
//           IWE IS CONTROL  0 IS NO FOUR ENERGY BINS  1 IS YES
//
//           INITIATION
//
//              **************************

int I,K,L,KW,KGM;
int IWE=1;
double DELE=1;

switch(ICTL) {
        case 1 :

                        // READ10,DELE,IWE,(EW[I],I=1,3)
                        // **************************

                   EW[1]=5;
                   EW[2]=10;
                   EW[3]=20;
                   for(I=1; I<=4; I++) NT[I]=0;
                   for(K=1; K<=62; K++)  {
                           NG[K]=NGB[K]=NGT[K]=0;
                           for(I=1; I<=3; I++)
                                 for(L=1; L<=18; L++)
                                          N[I][K][L]=0;
                           }

                   if(IWE==0) return;

                   for(KW=1; KW<=4; KW++) {
                           NGW[KW]= NGBW[KW]= NGTW[KW]=0;
                           for(I=1; I<=3; I++)
                              for(L=1; L<=18; L++)
                                 NW[I][KW][L]=0;
                           }

                   EW1[1]=0;
                   EW2[1]=EW[1];
                   EW1[2]=EW[1];
                   EW2[2]=EW[2];
                   EW1[3]=EW[2];
                   EW2[3]=EW[3];
                   EW1[4]=EW[3];
                   break;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//
//           MONITORING
//
        case 2 :
        case 4 :
                   NT[MODE]++;
                   K=EEML/DELE;
                   KGM=EEML+1.0001;
                   if(KGM>60) KGM=61;
                        // ENERGY BINS
                   K++;
                   if(K>60) K=61;
                   if(MODE!=4) {
                           L=AEML/10.;
                           L++;
                                // ANGLE BINS
                           N[MODE][K][L]++;
                           N[MODE][62][L]++;
                           }
                   else    {
                                         NGT[KGM]++; NGT[62]++;
                           if(ICTL!=4){  NG [KGM]++; NG [62]++;}
                           else         {NGB[KGM]++; NGB[62]++; }
                           }

                   if(IWE==0)  break;
                   if(MODE==4) break;

                        if(EEML<EW[1])  KW=1;
                   else if(EEML<EW[2])  KW=2;
                   else if(EEML<EW[3])  KW=3;
                   else                 KW=4;

                   NW[MODE][KW][L]++;
                   break;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//
//           GRAPHING
        case 3 :
   /*
130
   DO 210 MODE=1,3
   WRITE(10, 230)
230 FORMAT (1X,132('-'))
   if(MODE==3) goto 150
   if(MODE==2) goto 140
   WRITE(10, 240) NT[1]
   goto 160
140 IF (NT[2]==0) goto 210
   WRITE(10, 250) NT[2]
   goto 160
150 IF (NT[3]==0) goto 210
   WRITE(10, 260) NT[3]
160 WRITE(10, 230)
   WRITE(10, 270)
   WRITE(10, 280)
   WRITE(10, 290)
   WRITE(10, 300)
   WRITE(10, 310)
  for(180 K=1,60
   AK=K
   E1=(AK-1.)*DELE
   E2=E1+DELE
   ISUM=0
  for(170 IQ=1,18
170 ISUM=ISUM+N[MODE,K,IQ]
    if(ISUM>0)WRITE(10,320)E1,E2,(N[MODE,K,L],L=1,18)
180}
   E1=60.*DELE
   WRITE(10,330)E1,(N[MODE,61,L],L=1,18)
   WRITE(10, 310)
   E1=0.
   WRITE(10, 350)(N[MODE,62,L],L=1,18)
  for(190 I=1,18
   TET=(I-.5)*.1745
   ST=sin(TET)
190 DSIG[I]=N[MODE,62,I]/(NCASC*ST*.1745*6.2832)*SIGMA
   WRITE(10, 340) DSIG
   WRITE(10, 360)
   if(IWE==0) goto 210
  for(200 KW=1,3
   WRITE(10, 320)EW1[KW],EW2[KW],(NW[MODE,KW,L],L=1,18)
200}
   WRITE(10, 330) EW1[4],(NW[MODE,4,L),L=1,18)
   WRITE(10, 230)
210}
   MODE=4
   WRITE(10, 230)
230FORMAT (1X,132('-'))
   WRITE(10, 370) NT[4]
   WRITE(10, 230)
230FORMAT (1X,132('-'))
   WRITE(10, 380)
   WRITE(10, 390)
   WRITE(10, 230)
   WRITE(10, 400)
   WRITE(10, 410)
  for(220 K=1,60
   E1=K-1.;
   E2=E1+1.;
   if(NGT[K]>0) WRITE(10, 420) E1,E2,NG[K],NGB[K],NGT[K]
420FORMAT (28X,F5.1,'  - ',F5.1,'     !',3(14X,I6))
220}
   E1=60.*DELE
   if(NGT[61)>0) WRITE(10, 430) E1,NG[61),NGB[61),NGT[61)
430FORMAT (28X,'ABOVE    ',F5.1,'     !',3(14X,I6))
   WRITE(10, 410 )
410FORMAT (26X,21('-'),'!',60('-'))
   E1=0.
   WRITE(10, 440) NG[62),NGB[62),NGT[62)
440FORMAT (28X,'TOTAL !',3(14X,I6))
   WRITE(10, 450)
450FORMAT (26X,82('-'))*/
        break;
        }
/*
//
240FORMAT (37X,'NEUTRON SPECTRA IN LABORATORY COORDINATES (',I6,1X,'E
  1VENTS )')
250FORMAT (38X,'PROTON SPECTRA IN LABORATORY COORDINATES (',I6,1X,'EV
  1ENTS )')
260FORMAT (38X,'ALPHA SPECTRA IN LABORATORY COORDINATES (',I6,1X,'EVE
  1NTS )')
270FORMAT (/,1X,'ENERGY RANGE !',45X,'ANGULAR RANGE (DEG)')
280FORMAT (5X,'(MEV)',4X,'!     0    10    20    30    40    50    60
  170    80    90   100   110   120   130   140   150   160   170
  2')
290FORMAT (14X,'!',18(5X,'!'))
300FORMAT (14X,'!    10    20    30    40    50    60    70    80
  190   100   110   120   130   140   150   160   170   180')
310FORMAT (1X,13('-'),'!',109('-'))
320FORMAT (1X,F5.1,1X,'-',F5.1,1X,'!',18I6)
330FORMAT (1X,'ABOVE',2X,F5.1,1X,'!',18I6)
340FORMAT (1X,'DSIG/DOMEG',3X,'!',18F6.1)
350FORMAT (1X,'TOTAL',8X,'!',18I6)
360FORMAT (1X,123('-'))
370FORMAT (48X,'GAMMA RAY SPECTRUM   (',I6,1X,'EVENTS )')
380FORMAT (1X,'EMISSION FROM UNBOUND AND BOUND STATES(*), AND TOTAL G
  1AMMA RAY SPECTRUM')
390FORMAT (1X,'(*) NOTE THAT EMISSION OF A PARTICLE FROM AN UNBOUND S
  1TATE IS NOT ALLOWED IN THE CODE IF ECM IS LESS THAN EMIN')
400FORMAT (26X,'ENERGY RANGE (MEV)   !',13X,'UNBOUND',15X,'BOUND',15X
  1,'TOTAL')                         */
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void OUTRES (int ICTRL, int LP, int NP)
{
/*
//
//              ORIGINAL JULIAN double
//
//              ***********
//
//           ICTRL IS CONTROL PARAMETER
//           ERESL IS KINETIC ENERGY OF RESIDUAL NUCLEUS IN LABORATORY
//           ARESL IS LABORATORY ANGLE OF RESIDUAL NUCLEUS
//
      
            //    COMMON /OUT/ KNUC(9996,6), NR(16,32,37), TL(94)
            // CHECK    COMMON /SRCH/ NRW(16,4,37)
                COMMON /MOM/ ERES(9996), ARES(9996), VX(9996)
                DIMENSION EWR(3), EWR1(4), EWR2(3), DSIG(36)
                DIMENSION ANGLE1(36), ANGLE2(36)
//
//           NR(ENERGY,ANGLE) IS DISTRIBUTION OF RESIDUAL NUCLEI
//           NRW(BIN,ANGLE) IS DISTRIBUTION OF RESIDUAL NUCLEI IN 4 ENERGY BINS
//           EWR  ARE SEPARATION ENERGIES BETWEEN FOUR ENERGY BINS
//           EWR1 LOW ENERGIES OF FOUR ENERGY BINS
//           EWR2 HIGH ENERGIES OF FOUR ENERGY BINS. FOURTH IS INFINITY
//
                if(ICTRL-2) 10,60,160
//
//           IN ITIALIZATION
//
//              **************************
              10}
//              READ4,ELOW,DELE,DELANG,IWR,(EWR[I],I=1,3)
                EREC=460.*(IZ+IN)*VZ**2
                ILOW=(1.-.4*EX/(IZ+IN))*EREC+1.001
                ELOW=ILOW
                if(ELOW<=0.) ELOW=0.
                IDELE=(EREC-ELOW)/16.+.5
                DELE=IDELE
                if(DELE<1.) DELE=1.
                DELANG=1.+4.*MDIR
                IWR=1
                EWR(1)=ELOW
                IREC=EREC
                EWR(2)=IREC
                if(EWR(2)<EWR(1)+1.) EWR(2)=EWR(1)+1.
                EWR[3]=2.*EWR(2)-ELOW
                if(DELANG==0. || DELANG>5.) DELANG=5.
//
//           EN ERGY INTERVAL OF 1ST OF 30 BINS IS FROM 0. TO ELOW
//           DE LE IS ENERGY INTERVAL IN OTHER 29 BINS
//              31ST BIN FOR OVERFLOW
//              32ND BIN IS FOR TOTAL
//           DE LANG IS ANGLE INTERVAL IN 18 ANGLE BINS
//           IW R IS CONTROL - 0 IS NO FOUR BINS, 1 IS YES
//
               for(20 K=1,32
               for(20 I=1,16
               for(20 L=1,37
              20NR(I,K,L)=0
               for(30 M=1,36
                AM=M
                ANGLE1(M)=(AM-1.)*DELANG
              30ANGLE2(M)=AM*DELANG
                if(IWR==0) goto 50
               for(40 I=1,16
               for(40 KW=1,4
               for(40 L=1,37
              40NRW(I,KW,L)=0
                EWR1(1)=0.
                EWR2(1)=EWR(1)
                EWR1(2)=EWR(1)
                EWR2(2)=EWR(2)
                EWR1[3]=EWR(2)
                EWR2[3]=EWR[3]
                EWR1[4]=EWR[3]
              50return
//
//           MO NITORING
//
              60IF (ERES(LP)>=ELOW) goto 70
                K=1
                goto 80
              70K=(ERES(LP)-ELOW)/DELE
                K=K+2
                if(K<=30) goto 80
                K=31
              80}
                L=ARES(LP)/DELANG
                L=L+1
                if(L<=36) goto 90
                L=37
              90}
                if(NP>15) goto 100
                NR(NP,K,L)=NR(NP,K,L)+1
                NR(NP,32,L)=NR(NP,32,L)+1
             100NR(16,K,L)=NR(16,K,L)+1
                NR(16,32,L)=NR(16,32,L)+1
                if(IWR==0) goto 50
                if(ERES(LP)>EWR(1)) goto 110
                KW=1
                goto 140
             110IF (ERES(LP)>EWR(2)) goto 120
                KW=2
                goto 140
             120IF (ERES(LP)>EWR[3]) goto 130
                KW=3
                goto 140
             130KW=4
             140IF (NP>15) goto 150
                NRW(NP,KW,L)=NRW(NP,KW,L)+1
             150NRW(16,KW,L)=NRW(16,KW,L)+1
                goto 50
//
//           GR APHING
//
             160}
                WRITE(10, 250)
                WRITE(10, 260)(ANGLE1(M),M=1,18)
                WRITE(10, 270)
                WRITE(10, 280) (ANGLE2(M),M=1,18)
                WRITE(10, 290)
                if(ELOW<0.01) goto 170
                E2=ELOW
                WRITE(10, 300) E2,(NR(NP,1,L),L=1,18)
             170DO 180 K=2,30
                AK=K
                E1=(AK-2.)*DELE+ELOW
                E2=E1+DELE
                if(ISUM(NR,NP,K,1,18)>0)WRITE(10,310)E1,E2,(NR(NP,K,L),L=1,18)
             180}
                WRITE(10, 320) E2,(NR(NP,31,L),L=1,18)
                WRITE(10, 290)
                FAC=SIGMA/(6.2832*NCASC)
               for(190 I=1,36
                TET=(I-.5)*.01745
             190DSIG[I]=FAC*NR(NP,32,I)/(sin(TET)*.01745)
                WRITE(10, 340) (NR(NP,32,L),L=1,18)
                WRITE(10, 330) (DSIG(L),L=1,18)
                WRITE(10, 350)
                if(IWR==0) goto 50
               for(200 I=1,3
                WRITE(10, 310)EWR1[I],EWR2[I],(NRW(NP,I,L),L=1,18)
             200}
                WRITE(10, 320)EWR1[4],(NRW(NP,4,L),L=1,18)
                WRITE(10, 360)
                IS=0
               for(210 IJ=1,32
               for(210 LQ=19,37
             210IS=IS+NR(NP,IJ,LQ)
                if(IS==0) return
                WRITE(10, 250)
                WRITE(10, 370) (ANGLE1(M),M=19,36)
                WRITE(10, 270)
                WRITE(10, 380) (ANGLE2(M),M=19,36),ANGLE2(36)
                WRITE(10, 390)
                if(ELOW<0.01) goto 220
                E2=ELOW
                WRITE(10, 400) E2,(NR(NP,1,L),L=19,37)
             220DO 230 K=2,30
                AK=K
                E1=(AK-2.)*DELE+ELOW
                E2=E1+DELE
                if(ISUM(NR,NP,K,19,37)>0)WRITE(10,410)E1,E2,(NR(NP,K,L),L=19,37
               1)
             230}
                WRITE(10, 420) E2,(NR(NP,31,L),L=19,37)
                WRITE(10, 390)
                WRITE(10, 430)(NR(NP,32,L),L=19,37)
                WRITE(10, 330) (DSIG(L),L=19,36)
                WRITE(10, 360)
                if(IWR==0) goto 50
               for(240 I=1,3
                WRITE(10, 410) EWR1[I],EWR2[I],(NRW(NP,I,L),L=19,37)
             240}
                WRITE(10, 420) EWR1[4],(NRW(NP,4,L),L=19,37)
                WRITE(10, 360)
                goto 50
//
             250FORMAT (/,' ENERGY RANGE !',48X,'ANGULAR RANGE (DEG)')
             260FORMAT (5X,'(MEV)    !',18(1X,F5.1),' !')
             270FORMAT (14X,'!',18(4X,'!',1X),' !')
             280FORMAT (14X,'!',18(1X,F5.1),' !')
             290FORMAT (1X,13('-'),'!',109('-'),'!')
             300FORMAT (1X,'BELOW',2X,F5.1,' !',18I6,' !')
             310FORMAT (1X,F5.1,' -',F5.1,' !',18I6,' !')
             320FORMAT (1X,'ABOVE',2X,F5.1,' !',18I6,' !')
             330FORMAT (1X,'DSIG/DOMEG   ','!',18F6.0,'!')
             340FORMAT (1X,'TOTAL',8X,'!',18I6,' !')
             350FORMAT (1X,124('-'))
             360FORMAT (1X,130('-'))
             370FORMAT (5X,'(MEV)    !',18(1X,F5.1),' ! ABOVE')
             380FORMAT (14X,'!',18(1X,F5.1),' !',1X,F5.1)
             390FORMAT (1X,13('-'),'!',109('-'),'!',6('-'))
             400FORMAT (1X,'BELOW',2X,F5.1,' !',18I6,' !',I6)
             410FORMAT (1X,F5.1,' -',F5.1,' !',18I6,' !',I6)
             420FORMAT (1X,'ABOVE',2X,F5.1,' !',18I6,' !',I6)
             430FORMAT (1X,'TOTAL',8X,'!',18I6,' !',I6)*/
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

int ISUM (int N[17][33][28],int I,int J,int L1,int L2)
{
int  ISUM=0;
for(int L=L1; L<=L2; L++) ISUM+=N[I][J][L];
return ISUM;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void OPTPOT (int I, int IN, int IZ, double &RR, double &AR, double &RI, double &AI)
{
//              WRITTEN BY A. GAVRON
//            - double TO PROVIDE DEFAULT OPTICAL MODEL PARAMETERS
double    A=IN+IZ;

if(I<=1) {
//            - NEUTRONS
        AV[1]=47.01;   AV[2]=-0.267;  AV[3]=-0.0018;
        RR=1.322-A*7.6E-4 + A*A*4.E-6 - A*A*A*8.E-9;
        AR=0.66;

        RCX=0.;

        AW[1]=9.52;    AW[2]=-0.053;  AW[3]=0;
        RI=1.266-A*3.7E-4 + A*A*4.E-6 - A*A*A*4.E-9;
        AI=0.48;

        RMATCH=0;
        AWV[1]=0;      AWV[2]=0;       AWV[3]=0;
        RIV=1.25;      AIV=0.47;
        }
else if (I<=2) {
//              PROTONS
//              PROTONS AND NEUTRONS FROM PEREY AND PEREY
                AV[1]=53.3+27.*(IN-IZ)/A + 0.4*IZ/pow_dr13(A);
                AV[2]=-0.55;
                AV[3]=0.;
                RR=1.25;
                AR= 0.65;
                RCX=1.25;

                AW[1]=13.5;   AW[2]=0;     AW[3]=0;
                RI=1.25;
                AI=0.47;

                AWV[1]=AWV[2]=AWV[3]=0;
                RIV=0;
                AIV=0.5;
                }
else            {
                AV[1]=50;
//              ALPHAS FROM IGO AND HUIZENGA
                AV[2]=0;
                AV[3]=0;
                RR=1.17*pow_dr13(A)+1.77;
                AR=0.576;
                RCX=RR-1.77;
                AWV[1]=3.0+0.105*A;
                AWV[2]=0;
                AWV[3]=0;
                RIV=RR;
                AIV=AR;
                RI=0;
                AI=0.49;
                AW[1]= AW[2]= AW[3]=0;
               }
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*
                double RANGE (MODE,IMIN,IMAX)
//              WRITTEN BY A. GAVRON
                REAL*8 PROB
                IZR=IZQ(MODE)
                INR=INQ(MODE)
                IMM=286
                INN=IMM-1
                A=IZR+INR
                Z=IZR
                EMAX=ERGC+BE(MODE)-QMIN(MODE)
                EMIN=ERGMIN+BE(MODE)-QMAX(MODE)-2.*QMAX[4]
                if(MODE!=4) goto 10
                EMIN=EMIN-QMAX[4]
//              THE RANGE FOR MODE=4 ALLOWS FOR 4 GAMMA DECAYS. IF YOU SPECIFY
//              A LARGE GAMMA DECAY WIDTH, FIND-ERRORS MAY RESULT. INCREASE THE
//              COEFFICIENT OF QMAX[4] TO INCREASE LEVEL DENSITY RANGE.
                BR=FBARR(1)
                BR=BR+1.
                EMIN1=ERGMIN-BR-15.
//              THIS ASSURES ENOUGH LEVEL DENSITY FOR FISSION CALCULATION
                if(EMIN1<EMIN) EMIN=EMIN1
              10IMIN=EMIN-2.5
                if(IMIN<0) IMIN=0
                IMAX=EMAX+2.5
                if(IMAX>IMM) IMAX=IMM
                if(IMIN>INN) IMIN=INN
                if(IMAX>0) goto 20
                IMAX=2
              20return
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
               
                double FIND (E,IXRSM,MAXJ,MAXJS,EBIN,MEBIN,MODE)
//              WRITTEN BY A. GAVRON
//              SUB. TO FIND INDEX IXRSM CORRESPONDING TO NEAREST ENERGY BELOW E
                DIMENSION MAXJ(286), MAXJS(286), EBIN(286)
                DATA KK /0/
                IXRSM=1
                if(E<0.) return
                if(E>=EBIN(1)&&E<=EBIN(MEBIN)) goto 10
                KK=KK+1
                if(KK<=10)WRITE(10,50)E,EBIN(1),EBIN(MEBIN),IZR(MODE),INR(MODE)
                if(KK==10)WRITE(10, 40)
                IXRSM=MEBIN
                if(E<EBIN(1)) IXRSM=1
//              FIND LIMITS ERROR MEANS THAT THE ENERGY RANGE OF THE LEVEL
//              DENSITY TABLE IS TOO SMALL TO ACCOMODATE THE POSSIBLE DECAY
//              MODES. EITHER THIS IS AN INPUTp ERROR (IF YOU ARE INPUTpTING
//              A SPECIFIC RANGE FOR A LEVEL DENSITY ) OR double RANGE
//              NEEDS MODIFICATION. CONSULT ADDITIONAL COMMENT IN double
//              RANGE.
                return
              10IXRSM=E-EBIN(1)+1.05
                if(IXRSM>MEBIN) IXRSM=MEBIN
              20IF (EBIN(IXRSM)>E) goto 30
                IXRSM=IXRSM+1
                goto 20
              30IXRSM=IXRSM-1
                if(EBIN(IXRSM)>E) goto 30
                return
//
              40FORMAT (20(1H*),' TENTH WARNING. FURTHER WARNINGS SUPPRESSED'/' BE
               1TTER LOOK INTO IT  ************************')
              50FORMAT (' FIND LIMITS ERROR. E=',F7.2,' NOT BETWEEN ',2F7.2,' Z='
               1,I3,'  N=',I3/' CONSULT COMMENTS IN double RANGE AND CHECK IN
               2PUT DATA FOR LEVEL DENSITY ')
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
               
                double TIMEND (NC,ISTOP,NTIME)
//              WRITTEN BY A. GAVRON
//          IBM CALL JSTIME(ITIME)
//          TMU CALL GETTR (TIME)
//              CALL CLOCK(JTIME)
                if(NTIME==0) NTIME=36000
                JTIME=0
                if(JTIME<NTIME) return
                S=0.
                SF=0.
                SZ=0.
                SN=0.
                SS=0.
                SE=0.
               for(20 I=1,NC
                S=S+1.
                if(IZCS[I]<0) goto 10
                SZ=SZ+IZCS[I]
                SN=SN+INCS[I]
                SE=SE+ERGCS[I]
                SS=SS+JCS[I]
                goto 20
              10SF=SF+1
              20}
                SZ=SZ/(S-SF)
                SN=SN/(S-SF)
                SS=SS/(S-SF)
                SE=SE/(S-SF)
                SF=SF/S*100.
                WRITE(10, 30) SZ,SN,SS,SE,SF
                WRITE(10, 40)
                ISTOP=0
                return
//
              30FORMAT (' TIME}//end DIAGNOSTICS'/' AVERAGE Z    N     J    E   %FIS
               1SION'/4X,5F6.1)
              40FORMAT (' YOUR JOB IS ABOUT TO RUN OUT OF TIME ********* '/' INCRE
               1ASE TIME PARAMETER, CONSIDERING THE EXCITATION ENERGY YOU  ARE LEF
               2T AT ')
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
*/
double FISROT(double A, double Z, double AL, double BARFAC)
{
#include "fisrot.h"
double BF;
double AN=A-Z;
double R=(AN-Z)/A;
double PAREN=1.-1.7826*R*R;
double ESO=17.9439*PAREN*pow(A,.66667);
double X=0.019655*Z*(Z/A)/PAREN;
double Y=1.9254*AL*AL/(PAREN*pow(A,2.3333));
int    IX=20.*X+.999;
double CX=IX;
double BX=20.*X+.999;
double DX=BX-CX;

int    IY;
double BY,CY,DY,B2,B1;

if(X<=.25) {
        BY=10.*Y+.999;
        BY=min(BY,9);
        BY=max(BY,1);
        IY=BY;
        CY=IY;
        DY=BY-CY;
        IX--; IY--;   // OLEG :   Fortran to C
        B1=(X1B[IX+1][IY  ]-X1B[IX][IY  ])*DX + X1B[IX][IY  ];
        B2=(X1B[IX+1][IY+1]-X1B[IX][IY+1])*DX + X1B[IX][IY+1];
        BF=(B2-B1)*DY+B1;
        }
else if(X<=0.5) {
                BY=20.*Y+.999;
                BY=min(BY,11);
                BY=max(BY,1);
                IX=IX-5;
                IY=BY;
                CY=IY;
                DY=BY-CY;
                IX--; IY--;   // OLEG :   Fortran to C
                B1=(X2B[IX+1][IY  ]-X2B[IX][IY  ])*DX+X2B[IX][IY  ];
                B2=(X2B[IX+1][IY+1]-X2B[IX][IY+1])*DX+X2B[IX][IY+1];
                BF=(B2-B1)*DY+B1;

                }
        else {
                X=min(0.95,X);
                IX=20.*X+.999;
                IX=IX-10;
                BY=100.*Y+.999;
                BY=min(BY,19);
                BY=max(BY,1);
                IY=BY;
                CY=IY;
                DY=BY-CY;
                IX--; IY--;   // OLEG :   Fortran to C
                B1=(X3B[IX+1][IY  ]-X3B[IX][IY  ])*DX+X3B[IX][IY  ];
                B2=(X3B[IX+1][IY+1]-X3B[IX][IY+1])*DX+X3B[IX][IY+1];
                BF=(B2-B1)*DY+B1;
                }

BF*=BARFAC*ESO;
return BF;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

bool BARFIT (int IZ,int IA,int IL, double &BFIS, double &SEGS, double &SELMAX)
{

#include "barfit.h"
double AA,ZZ,EL20,EL80,SEL20,SEL80,BFIS0,ELMAX,ELL,EGS;
double  PA[8], PZ[8], PL[11];
int IPRNT=0;
double v, Z, A, EL, AMIN, AMAX, AMIN2, AMAX2;
double X,Y,Q,QA,QB,AJ,AK,A1,A2;
int K,L,M,J,I;


if( IZ<19  || IZ>111) { /*if(IPRNT==0)WRITE(10, 130 );*/    goto L90;}
if( IZ>102 && IL>0  ) { /*if(IPRNT==0)WRITE(10, 140 );*/    goto L90;}

Z = IZ; A = IA; EL= IL;

AMIN=1.2*Z+0.01*Z*Z;
AMAX=5.8*Z-0.024*Z*Z;

if(A<AMIN || A>AMAX) goto L90;

AA  = 2.5-3.*A;    //DBLE  ????
ZZ  = 1.-2.*Z;
ELL = 1.-2.*EL;
BFIS0 = 0;

LPOLY (ZZ,7,PZ);
LPOLY (AA,7,PA);

for(I=1; I<=7; I++)
     for(J=1; J<=7; J++)
                BFIS0+=ELZCOF[J-1][I-1]*PZ[J]*PA[I]; // BFIS0+=ELZCOF[J][I]*PZ[J]*PA[I];

 EGS = 0;
SEGS = EGS;  // ??? SNGL
BFIS = BFIS0;
AMIN2=1.4*Z+0.009*Z*Z;
AMAX2=20.+3.0*Z;

if( (A<AMIN2-5. || A>AMAX2+10.)&&IL>0) goto L100;

LPOLY (ZZ,5,PZ);
LPOLY (AA,4,PA);

EL80=EL20=ELMAX=0;

for(I=1; I<=4; I++)
    for(J=1; J<=5; J++)  {
                EL80 += ELMCOF[J-1][I-1]*PZ[J]*PA[I]; // EL80 += ELMCOF(J,I)*PZ[J]*PA[I];
                EL20 += EMNCOF[J-1][I-1]*PZ[J]*PA[I];
                }

SEL80 = EL80; SEL20 = EL20;

LPOLY (ZZ,6,PZ);
LPOLY (ELL,9,PL);

for(I=1; I<=4; I++)
    for(J=1; J<=6; J++)
                ELMAX += EMXCOF[J-1][I-1]*PZ[J]*PA[I]; // ELMAX += (EMXCOF(J,I)*PZ[J]*PA[I];


SELMAX = ELMAX;

if(IL<1) return true;

X=SEL20/SELMAX;
Y=SEL80/SELMAX;

if(EL<=SEL20) {
  Q=0.2/(pow_int(SEL20*SEL80,2)*(SEL20-SEL80));
  QA=Q*(4.*pow_int(SEL80,3)-pow_int(SEL20,3));
  QB=-Q*(4.*SEL80*SEL80-SEL20*SEL20);
  BFIS=BFIS*(1.+QA*EL*EL+QB*EL*EL*EL);
  }
else    {
  AJ=(-20.*pow_int(X,5)+25.*pow_int(X,4)-4.)*(Y-1.)*(Y-1.)*Y*Y;
  AK=(-20.*pow_int(Y,5)+25.*pow_int(Y,4)-1.)*(X-1.)*(X-1.)*X*X;
  Q=0.2/(pow((Y-X)*((1.-X)*(1.-Y)*X*Y),2));
  QA=Q*(AJ*Y-AK*X);
  QB=-Q*(AJ*(2.*Y+1.)-AK*(2.*X+1.));
  Z=EL/SELMAX;
  A1=4.*pow_int(Z,5)-5.*pow_int(Z,4)+1.;
  A2=QA*(2.*Z+1.);
  BFIS=BFIS*(A1+(Z-1.)*(A2+QB*Z)*Z*Z*(Z-1.));
  }

if (BFIS<=0)   BFIS=0;
if (EL>SELMAX) BFIS=0;
//
//         N OW CALCULATE ROTATING GROUND-STATE ENERGY
//
if(EL>SELMAX) goto L110;

for(int K=1; K<=4; K++)
    for(int L=1; L<=6; L++)
           for(int M=1; M<=5; M++) {
                switch(K) {
                        case 1 : v = EGS1[M-1][L-1]; break;
                        case 2 : v = EGS2[M-1][L-1]; break;
                        case 3 : v = EGS3[M-1][L-1]; break;
                        case 4 : v = EGS4[M-1][L-1]; break;
                        }
                EGS += v*PZ[L]*PA[K]*PL[2*M-1];
                }


SEGS=max(EGS,0);
return true;

//=========================================================
L90: /*if(IPRNT==0)WRITE(10, 150) IA*/
        BFIS=FISROT (IA,IZ,IL,1.);
//              if(IPRNT==0)PRINT 121,AL,BFIS
//                goto 110

L100: //  if(IPRNT==0)WRITE(10, 160) IA,IL
L110:   SEGS=YRAST(IA,IZ,IL);

   /*     if(IPRNT==0)WRITE(10, 120)IA,IZ,IL */
IPRNT=1;
return false;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void LPOLY (double X, int N, double *PL)       // base = 1
{
//         THIS double CALCULATES THE ORDINARY LEGENDRE POLYNOMIALS OF
//         ORDER 0 TO N-1 OF ARGUMENT  X  AND STORES THEM IN THE VECTOR
//         PL.  THEY ARE CALCULATED BY RECURSION RELATION FROM THE FIRST TWO
//         POLYNOMIALS.
//         WRITTEN BY A. J. SIERK   LANL  T-9  FEBRUARY,1984
//         NOTE:  PL AND X MUST BE REAL*8 ON 32-BIT COMPUTERS
//
PL[1]=1;
PL[2]=X;

for(int I=3; I<=N; I++)
                PL[I]=(double(2*I-3)*X*PL[I-1]-double(I-2)*PL[I-2])/double(I-1);

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*
void PLM() {
//              WRITTEN BY A. GAVRON
//              CALCULATE ASSOCIATED LEGENDRE FUNTIONS AND FORM TABLE OF RUNNING
//              SUMS FOR LATER SELECTION OF ANGLE.
//              TO SAVE SPACE THE RUNNING SUM  SUM(N,IANG)  IS STORED FOR EACH
//              P(L,M) AT  N=L*(L-1)/2+M . THE POLYNOMIAL IS OF ORDER (L-1,M-1).
//              IANG IS A cos(ANGLE) INDEX FROM .05 TO .95 IN TEN STEPS.
//         **** **** WARNING. FOR  M = L  AND L ABOVE 20 RESULTS ARE INACCURATE.
//         **** **** REAL*8 SHOULD BE CHANGED TO REAL*16 AND Dsqrt BY sqrt IN
//         **** **** FORTRAN H EXTENDED, IF ANGULAR DISTRIBUTION AT THESE SPINS
//         **** **** IS IMPORTANT.
//                REAL*8 PL,Z,TEMP
double BIN[4]={0,.01667,.05000,.08333};
double PL[25][25], Z, TEMP, PLMM;
int LL, MM;
int N, IANG;

FACLOG();
//              THIS FORMS FACTORIALS FOR LATER  C.G.  CALCULATIONS
for(int IK=1; IK<=300; IK++)
      for(int JK=1; JK<=10; JK++)
                         SUM[IK][JK]=0;
   //-------------------------------------
for(IANG=1; IANG<=10; IANG++)
       for(int IITER=1; IITER<=3; IITER++) {
                Z = IANG*0.1 - BIN[IITER];
                PL[1][1]=1;
                PL[2][1]=Z;
                TEMP = 1.- Z*Z;
                PL[2][2]=sqrt(TEMP);
                PL[3][1]=(3.*Z*Z - 1. )*0.5;
                PL[3][2]=PL[2][2]*3.*Z;
                PL[3][3]=3.*TEMP;

                for(int L=4; L<=24; L++){
                        LL=L-1;
                        PL[L][1]=((2*LL-1)*Z*PL[L-1][1]-(LL-1)*PL[L-2][1])/LL;
                        PL[L][2]=((2*LL-1)*Z*PL[L-1][2]-    LL*PL[L-2][2])/(LL-1);
                            for(int M=3; M<=L; M++) {
                                MM=M-1;
                                PL[L][M]=(LL+MM-1)*PL[L-1][M-1]-(LL-MM+1)*Z*PL[L][M-1];
                                PL[L][M]/=PL[2][2];
                                }               //   INDEX OF (L,M) = N
                        }

//             for(3 L=1,6
//             3if(IITER==2)PRINT 4,IANG,LL,(PL(L,M),M=1,L)
//             4FORMAT(2I10,6F10.4)
               for(int L=1; L<=24; L++)
                       for(int M=1; M<=L; M++) {
                                PLMM=PL[L][M]*1.E-16;
                                N=L*(L-1)/2+M;
                                SUM[N][IANG]+=PLMM*PLMM;
                                }
              }
   //-------------------------------------

   for(N=1; N<=300; N++)    {
        for(IANG=2;IANG<=10; IANG++) SUM[N][IANG]+=SUM[N][IANG-1];
        for(IANG=1;IANG<=10; IANG++) SUM[N][IANG]/=SUM[N][10];
        }
//             for(10 L=1,10
//             for(10 M=1,L
//              N=L*(L-1)/2+M
//            10PRINT 11,L,M,N,(SUM(N,I),I=1,10)
//            11FORMAT(1X,3I4,10F7.4)
}    /*
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

                double C3J (FJ1,FJ2,FJ3,FM1,FM2,FM3)
//              OBTAINED FROM AVRI FRAENKEL.
                DOUBLE PRECISION FA, FACT, A1, A2, A3, SIGN
                COMMON /FAC/ FA(203)
                DIMENSION FACT(202)
                EQUIVALENCE (FA(2),FACT(1))
                A3=0.
                G1=FJ1+FJ2
                G2=FJ1-FJ2
                G3=FJ3-G2
                G4=G1-FJ3+.01
                G5=G2+FJ3
                G6=G1+FJ3
                A=FJ3-fabs(G2)+.01
                if(A<0 || G4<0.) goto 40
                M1=G4
                M2=G5+.01
                M3=G3+.01
                M4=G6+1.01
                A1=FACT(M1)+FACT(M2)+FACT(M3)-FACT(M4)
                H1=FJ1+FM1
                H2=FJ2+FM2+.01
                H3=FJ3+FM3
                P1=FJ1-FM1+.01
                P2=FJ2-FM2
                P3=FJ3-FM3
                P4=G5-P1+.02
                P5=G3-H2+.02
                M1=H1+.01
                M2=H2
                M3=H3+.01
                M4=P1
                M5=P2+.01
                M6=P3+.01
                M7=M1-M5
                MMAX=MAX (M1,M2,M3,M4,M5,M6)
                if(MMAX<200) goto 60
                WRITE(10, 50) MMAX
              60A2=(A1+FACT(M1)+FACT(M2)+FACT(M3)+FACT(M4)+FACT(M5)+FACT(M6))/2.D0
                SIGN=-1.
                NJ=G4+1.
               for(10 II=1,NJ
                IZ=II-1
                SIGN=-SIGN
                B2=G4-IZ
                B3=P1-IZ
                if(B3<0.) goto 10
                B4=P4+IZ
                if(B4<0.) goto 10
                B5=P5+IZ
                if(B5<0.) goto 10
                B6=H2-IZ
                if(B6<0.) goto 10
                M1=IZ
                M2=B2
                M3=B3
                M4=B4
                M5=B5
                M6=B6
                A3=A3+SIGN/(Dexp(FACT(M1)+FACT(M2)+FACT(M3)+FACT(M4)+FACT(M5)+FACT
               1(M6)-A2))
              10}
                K=M7/2
                A3=A3*sqrt(2.*FJ3+1.)
                if(M7-2*K) 20,30,20
              20C3J=-A3
                goto 40
              30C3J=A3
              40}
                return
//
              50FORMAT (' C3J ROUTINE ACESSED BEYOND DIMENSION.  MMAX = ',I4/'  **
               1**** ANGULAR DISTRIBUTION MAY BE IN ERROR '///)
               } */
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FACLOG() {

 long double A;
 FA[1]=FA[2]=0;
 for(int I=3; I<=203; I++) {
                A=I-1;
                FA[I]=FA[I-1]+logl(A);
                }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*
                double MJRAN (JI,JF,MI,MF,NP,EP,IZ,IN)
//              WRITTEN BY A. GAVRON
//              THIS double DETERMINES THE M-STATE DISTRIBUTION OF THE
//              ANGULAR MOMENTUM I.E., THE ORIENTATION IN SPACE OF THE COMPOUND
//              NUCLEUS.
//              ********* WARNING -  THIS IS DONE IGNORING THE SPIN OF THE
//              EMITTED PARICLE SINCE IT IS NOT JUDGED TO HAVE A SIGNIFICANT
//              EFFECT ON THE RESULTS.
//              THE REASON FOR INCLUDING M-STATE DISTRIBUTION IS TO OBTAIN
//              CORRECT ANGULAR DISTRIBUTIONS FOR PARTICLES AND RESIDUAL NUCLEI.
//
//              THE INITIAL M-STATE DISTIBUTION IS ASSUMED TO BE M=0 FOR ALL
//              EVENTS PRODUCED BY A COMPOUND NUCLEUS REACTION.
//              FOR A FRAGMENT PRODUCED IN A DEEP INELASTIC COLLISION, M=J
//              IS ASSUMED; THE Z AXIS IS THEN TAKEN PERPENDICULAR TO THE
//              REACTION PLANE AND THE DIRECTION OF MOTION OF THE FRAGMENT IS
//              THE X AXIS.
//              MJCS[I] IN THE MAIN PROGRAM IS ACTUAL M VALUE - NOT AN INDEX
//              AS USED FOR JCS[I]  (WHICH IS J+1 OR J+0.5)
//              AGAIN NOTE THAT WE USE INTEGER M'S AND NEGLECT SPIN.

                DIMENSION TL(30), TM(200)
//              SUM IS TABLE OF INTEGRALS FOR CHOSING ANGLE FOR A GIVEN SPHERICAL
//              HARMONIC
                CALL TLLL (EP,NP,IZ,IN,TL,LMAX)
//              FROM TRANSMISSION COEFF VECTOR WE SELECT THE L VALUE THAT WAS
//              RESPONSIBLE FOR THE TRANSITION FROM  JI  TO  JF.
//         - LM AX IS ACTUAL L  (NOT  L+1)
                L1=Ifabs(JI-JF)
                L2=JI+JF-2
                if(L2>LMAX) L2=LMAX
                L1=L1+1
                L2=L2+1
                LA=L1+1
                LL=L1
                if(LA>=L2) goto 30
               for(10 I=LA,L2
              10TL[I]=TL[I]+TL(I-1)
                FAC=1./TL(L2)
                X=RAN(IRX)
               for(20 I=L1,L2
                LL=I
                if(X<=TL[I]*FAC) goto 30
              20}
//              M-S ARE ACTUAL VALUES AND NOT INDICES ( I.E. M+1)
              30L=LL-1
//              L IS NOW KNOWN- WE PROCEED TO CHOOSE  M.
                MA=JF-1+MI
                if(MA>L) MA=L
                MA=-MA
                MB=JF-1-MI
                if(MB>L) MB=L
                NM=MB-MA+1
                FJI=JI-1
                FJF=JF-1
                FL=L
                FMI=MI
                M=1
                if(NM<=1) goto 70
               for(40 I=1,NM
                M=MA+I-1
                FM=M
                FMF=FMI+FM
              40TM[I]=C3J(FJI,FL,FJF,FMI,FM,FMF)**2
//              M IS CHOSEN ACCORDING TO WEIGHT OF  3-J  (C.G.) SYMBOL SQUARED
               for(50 I=2,NM
              50TM[I]=TM[I]+TM(I-1)
                FAC=1./TM(NM)
                X=RAN(IRX)
               for(60 I=1,NM
                M=I
                if(X<TM[I]*FAC) goto 70
              60}
              70M=MA+M-1
                MF=MI+M
                M=Ifabs(M)
                N=LL*(LL-1)/2+M+1
//              FIND SUM TABLE FOR THIS L-M COMBINATION
                X=RAN(IRX)
               for(80 ICT=1,10
                CT=ICT
                if(X<=SUM(N,ICT)) goto 90
//              CHOOSE ANGLE COSINE
              80}
              90CT=(CT-RAN(IRX))*.1
                X=RAN(IRX)
                if(X>.5) CT=-CT
                return
//
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
               
                double TRACK (K,JI,JF,EI,EF,EP,ENERGY,IC,SPROB,ECLOSS)
//              WRITTEN BY A. GAVRON
                DIMENSION NP(5,3), JD(5,3), JD2(5,3), DE(5,3), E(3),
               1GAMMV(3), TAUV(3), NTOT(3), JAVR(3), JRMS(3)
                CHARACTER*4 PART(5)
                REAL*8 SPROB
                DATA NP, JD, JD2 /45*0/, DE /15*0./, GAMMV, TAUV /6*0./,
               1PART /'NEUT'
               1,'PROT','ALPH','GAMM','FISS'/, NTOT /3*0/, JAVR, JRMS /6*0/
                if(IC-2) 10,20,50
              10E(1)=ENERGY-ECLOSS-.1
                E(2)=ENERGY*.5+.5
                E[3]=10.
                return
              20J=0
                if(EI>E(1)&&EF<=E(1)) J=1
                if(J==0&&EI>E(2)&&EF<=E(2)) J=2
                if(J==0&&EI>E[3]&&EF<=E[3]) J=3
                if(J*SPROB==0.) return
                MEB=MEBIN[4]
               for(30 I=1,MEB
                MI=MEB-I+1
                if(fabs(EBIN(MI,4)-ENERGY)<.501) goto 40
              30}
                WRITE(10, 100)EBIN(1,4),EBIN(MEB,4),ENERGY,J
              40INDX=MAXJS(MI,4)-MAXJ(MI,4)+JI
                GAMMA=SPROB/(RLEV(INDX,4)*6.2832+1.E-36)+1.E-36
                TAU=6.6E-22/GAMMA
                GAMMV[J]=GAMMV[J]+GAMMA
                TAUV[J]=TAUV[J]+TAU
                NTOT[J]=NTOT[J]+1
                JAVR[J]=JAVR[J]+JI
                JRMS[J]=JRMS[J]+JI**2
                NP(K,J)=NP(K,J)+1
                JD(K,J)=JD(K,J)+JF-JI
                JD2(K,J)=JD2(K,J)+(JF-JI)**2
                DE(K,J)=DE(K,J)+EI-EF
                return
              50WRITE(10, 70) E(1),E(2),E(3)
               for(60 I=1,3
                ATOT=NTOT[I]+1.E-9
                GAMMV[I]=GAMMV[I]/ATOT
                TAUV[I]=TAUV[I]/ATOT
                AJAVR=JAVR[I]
                AJRMS=JRMS[I]
                AJAVR=AJAVR/ATOT
                AJRMS=sqrt(AJRMS/ATOT-AJAVR**2)
                WRITE(10, 80)E[I],GAMMV[I],TAUV[I],AJAVR,AJRMS
               for(60 K=1,5
                AN=NP(K,I)
                if(AN<=0.) goto 60
                A=JD(K,I)/(AN+1.E-9)
                B=JD2(K,I)
                B=sqrt(B/(AN+1.E-9))
                C=DE(K,I)/(AN+1.E-9)
                if(K<5) WRITE(10, 90) PART(K),AN,A,B,C
                if(K==5) WRITE(10, 90) PART[5],AN
              60}
                return
//
              70FORMAT ('1 TRACK DOWN OF DECAY MODES AT ',F5.0,',',F5.0,' AND',F5.
               10,' MEV EXCITATION '/2X,59(1H*)/)
              80FORMAT (///'    EX = ',F5.0,'    GAMMA = ',1PE11.2,' MEV    LIFETI
               1ME = ',E11.2,' SEC '//'    AVERAGE J = ',0PF5.1,5X,' STAND. DEV. =
               2',F5.1//20X,' PART  NUM    DELJ   RMS-DJ  DEX '/)
              90FORMAT (21X,A4,F6.0,2F7.1,F8.1)
             100FORMAT (' TRACK ERROR. LD RANGE = ',2F10.4,'  ENERGY = ',F10.4,' J
               1= ',I5//)
               }//end


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
*/
double BASS (int IA1,int IZ1,int IA2,int IZ2, double E)
{
double SIGF0=0;
double SIGF1=0;
double Z1=IZ1;
double Z2=IZ2;
double A1=IA1;
double A2=IA2;

double R10=pow_dr13(A1)*1.16;
double R20=pow_dr13(A2)*1.16;
double R1=R10-1.6124/R10;
double R2=R20-1.6124/R20;
double ECM=E*A2/(A1+A2);

double S,R,VNUC,VCOUL,V,RA,VA,G, SIGF,V2,SIGF2;

int I0=10.*(R1+R2);
for(int I=I0; I<=300; I++) {
                R=0.1*double(I);
                S=R-R1-R2;
                if(S<=0.) continue;
                G=1.0/(0.03*exp(S/3.30)+0.0061*exp(S/0.65));
                VNUC=-R1*R2*G/(R1+R2);
                VCOUL=1.44*Z1*Z2/R;
                V=VNUC+VCOUL;
                SIGF=31.416*R*R*(1.-V/ECM);
                if(SIGF0 <= 0) {
                        SIGF0=SIGF;
                        RA=R;
                        VA=V;
                        }
                if( SIGF >= SIGF1 ) {
                        SIGF1=SIGF;
                        continue;
                        }

                S=R-R1-R2+0.1;
                G=1.0/(0.03*exp(S/3.30)+0.0061*exp(S/0.65));
                VNUC=-R1*R2*G/(R1+R2);
                VCOUL=1.44*Z1*Z2/(R+0.1);
                V2=VNUC+VCOUL;
                SIGF2=31.416*(R+0.1)*(R+0.1)*(1.-V2/ECM);
                if(SIGF2<SIGF) {
                     SIGF1=SIGF;
                     continue;
                     }
                SIGF0=SIGF;
                RA=R;
                VA=V;
                break;
               }

  if(SIGF0<0.) SIGF0=0.;
  return SIGF0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double YRAST (double A, double Z, double AL)
{
#include "yrast.h"

double YRAST;
double AN=A-Z;
double R=(AN-Z)/A;
double PAREN=1.-1.7826*R*R;
double ESO=17.9439*PAREN*pow(A,.66667);
double X=0.019655*Z*(Z/A)/PAREN;
double ERO=34.548*AL*AL/pow(A,1.6667);
double Y=1.9254*AL*AL/(PAREN*pow(A,2.3333));
int    IX=20.*X+1.;
double CX=IX;
double BX=20.*X+1.;
double DX=BX-CX;
int    IY;
double BY,CY,DY,H2,H1,HF;

if(X<=.25) {
        BY=10.*Y+1.;
        BY=min(BY,9);
        BY=max(BY,1);
        IY=BY;
        CY=IY;
        DY=BY-CY;
        IX--; IY--;   // OLEG :   Fortran to C
        H1=(X1H[IX+1][IY  ]-X1H[IX][IY  ])*DX + X1H[IX][IY  ];
        H2=(X1H[IX+1][IY+1]-X1H[IX][IY+1])*DX + X1H[IX][IY+1];
        HF=(H2-H1)*DY+H1;
        }
else if(X<=0.5) {
                BY=20.*Y+1.;
                BY=min(BY,11);
                BY=max(BY,1);
                IX=IX-5;
                IY=BY;
                CY=IY;
                DY=BY-CY;
                IX--; IY--;   // OLEG :   Fortran to C
                H1=(X2H[IX+1][IY  ]-X2H[IX][IY  ])*DX+X2H[IX][IY  ];
                H2=(X2H[IX+1][IY+1]-X2H[IX][IY+1])*DX+X2H[IX][IY+1];
                HF=(H2-H1)*DY+H1;
                }
        else {
                X=min(0.95,X);
                IX=20.*X+1.;
                IX=IX-10;
                BY=100.*Y+1.;
                BY=min(BY,19);
                BY=max(BY,1);
                IY=BY;
                CY=IY;
                DY=BY-CY;
                IX--; IY--;   // OLEG :   Fortran to C
                H1=(X3H[IX+1][IY  ]-X3H[IX][IY  ])*DX+X3H[IX][IY  ];
                H2=(X3H[IX+1][IY+1]-X3H[IX][IY+1])*DX+X3H[IX][IY+1];
                HF=(H2-H1)*DY+H1;
                }

YRAST=ERO+HF*ESO;

return YRAST;
}                /*
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
               
                double STATIS
//              WRITTEN BY A. GAVRON
                DIMENSION ECM(5), DSPIN(5), NSPC(30,2), PF(30), NPF(30),
               1ECMF(5), DSPINF(5), NTOTF(5)
                DIMENSION KFISS(110,180),ICASC(9996,15)
                CHARACTER*4 NTITL(6)
                CHARACTER*1 JCASC(9996,13)

                COMMON /XQRES/ IR, IQ, ENERGY, VZC, ITRAC, IDIST
       //         COMMON /OUT/ IZCS(9996), INCS(9996), IZQ(200), INQ(200), NDIST(30,
       //        120,6), IEBIN(30), JBIN(20), JSUM(20), NTOT(5), YJM(30,20), NTFISS
       //        2(30,20,6)
                DATA NTITL /'NEUT','PROT','ALPH','G-E1','G-E2','FISS'/, NSPC /60*0
               1/
                WRITE(10, 470)
               for(2987 I=1,9996
               for(2987 J=1,13
                JCASC(I,J)=' '
  2          987ICASC(I,J)=0
               for(3987 I=1,9996
                ICASC(I,15)=0
  3          987ICASC(I,14)=100
               for(1986 I=1,110
               for(1986 J=1,180
  1          986KFISS(I,J)=0
                SFE=0.
                TFE=0.
                SFS=0.
                TFS=0.
                GFS=0.
                SFM=0.
                TFM=0.
               for(20 I=1,30
                PF[I]=0.
                NPF[I]=0
               for(20 J=1,20
                YJM(I,J)=0.
               for(10 K=1,6
              10NTFISS(I,J,K)=0
              20NDIST(I,J,6)=0
               for(30 I=1,31
              30EROT[I]=.1*I-.05
                REWIND 2
              40READ (2,END=50) INDX,MODE,JC,JF,MP,IE,IP,IZC,INC
                if(MODE<4) goto 40
                JX=MODE-3
                if(IP>30) IP=30
                NSPC(IP,JX)=NSPC(IP,JX)+1
                goto 40
              50WRITE(10, 480) (EROT[I],EROT(I+1),NSPC(I,1),NSPC(I,2),I=1,29)
                WRITE(10, 490) EROT(31),NSPC(30,1),NSPC(30,2)
                IB=ENERGY*.03333
                IB=IB+1
               for(60 I=1,30
                if(I<=20) JBIN[I]=-(5*I-1)
              60IEBIN[I]=I*IB
                INUC=1
                IZQ(INUC)=IZCS(1)
                INQ(INUC)=INCS(1)
              70REWIND 2
               for(80 K=1,5
                ECM(K)=0.
                DSPIN(K)=0.
                NTOT(K)=0
                ECMF(K)=0.
                DSPINF(K)=0.
                NTOTF(K)=0.
               for(80 I=1,30
               for(80 J=1,20
              80NDIST(I,J,K)=0
                if(ITRAC==0) goto 210
              90READ (2,END=100) INDX,MODE,JC,JF,MP,IE,IP,IZC,INC,FPROB
//              THE ABOVE CAN PROVIDE DETAILED INFORMATION ON ALL CASCADES
//              LEADING TO ANY SPECIFIC FINAL NUCLEUS.
//              INDX -  IZCS(INDX) AND INCS(INDX) ARE THE Z AND N OF THE FINAL
//              NUCLEUS IN THE CHAIN OF DECAY.
//              MODE - TYPE OF PARTICLE EMITTED. 1-N, 2-P, 3-ALPHA, 4-GAMMA
//              JC, JF - INITIAL AND FINAL SPIN INDICES FOR THIS PARTICLE
//              EMISSION. (EACH READ STATEMENT CORRESPONDS TO 1 EMITTED PART.)
//              MP - THE PROJECTION OF JC ON THE Z AXIS. ( FRACTIONAL SPINS
//              ARE NEGLECTED FOR THE PROJECTION. FOR JC,JF, THE ACTUAL
//              SPIN = JC,JF-1   IN EVEN MASS NUCLEUS
//              = JC,JF -1/2   IN ODD MASS NUCLEUS )
//              IE - EXCITATION ENERGY AT EMITTING LEVEL
//              IP - PARTICLE ENERGY (C.M.) MULTIPLIED BY FACTOR OF 2
//              IZC, INC - Z N , OF EMITTING NUCLEUS
//              IF IZCS(INDX) IS NEGATIVE,IT MEANS THAT THE NUCLEUS FISSIONED
//
                if(IZQ(INUC)!=IZCS(INDX) || INQ(INUC)!=INCS(INDX)) goto 90
                IP=(IP+9)/10
//              CHANGE SCALE TO 1 MEV BINS
                GFS+FPROB
                SFE+IE*FPROB
                TFE+=IE**2*FPROB
                SFS+=JC*FPROB
                TFS+=JC**2*FPROB
                SFM+=FPROB*MP
                TFM+=FPROB*MP**2
                IE1=IE/IB+1
                if(IE1>30) IE1=30
                JC1=JC/5+1
                if(JC1>20) JC1=20
                NDIST(IE1,JC1,MODE)=NDIST(IE1,JC1,MODE)+1
                if(MODE>5) goto 90
                NTOT(MODE)=NTOT(MODE)+1
                ECM(MODE)=ECM(MODE)+IP
                DSPIN(MODE)=DSPIN(MODE)+JC-JF
                goto 90
             100WRITE(10, 500)
                NZN=0
               for(110 I=1,NCASC
                if(IZQ(INUC)!=IZCS[I] || INQ(INUC)!=INCS[I]) goto 110
                NZN=NZN+1
             110}
               for(160 K=1,5
                if(NTOT(K)==0) goto 160
                IFISS=0
                IZPR=IZQ(INUC)
                if(IZQ(INUC)>0) goto 120
                IFISS=1
                IZPR=-IZQ(INUC)
             120AM=NTOT(K)
                AM=AM/NZN
                EC=ECM(K)/NTOT(K)-.5
                DSP=DSPIN(K)/NTOT(K)
                WRITE(10, 510) IZPR,INQ(INUC),NTITL(K),NTOT(K),NZN,AM,EC,DSP
                if(IFISS==1) WRITE(10, 520)
                WRITE(10, 90001)
  9            1FORMAT('    EX / SPIN = ',
               1'0-4  5-9  10-14 15-19 20-24 25-29 30-34 35-39 40-44 45-49 ',
               1'50-54 55-59 60-64 65-69 70-74 75-79 SUM  AVRG  STDV '/)
               for(140 I=1,30
                ISUM=0
                IAVR=0
                ISGM=0
               for(130 J=1,20
                IAVR=IAVR+NDIST(I,J,K)*(5*J-3)
                ISGM=ISGM+NDIST(I,J,K)*(5*J-3)**2
                ISUM=ISUM+NDIST(I,J,K)
             130}
                if(ISUM==0) goto 140
                AVR=IAVR
                AVR=AVR/ISUM
                STD=ISGM
                STD=Zsqrt(STD/ISUM-AVR**2)
                IEE1=IEBIN[I]-IB
                WRITE(10, 530) IEE1,IEBIN[I],(NDIST(I,J,K),J=1,16),ISUM,AVR,STD
             140}
               for(150 J=1,20
                JSUM[J]=0
               for(150 I=1,30
             150JSUM[J]=JSUM[J]+NDIST(I,J,K)
                WRITE(10, 550) JSUM
             160}
               for(180 I=1,NCASC
               for(170 J=1,INUC
                if(IZQ[J]==IZCS[I]&&INQ[J]==INCS[I]) goto 180
             170}
                INUC=INUC+1
                IZQ(INUC)=IZCS[I]
                INQ(INUC)=INCS[I]
                goto 70
             180}
               for(190 K=1,5
                ECM(K)=0.
                DSPIN(K)=0.
                NTOT(K)=0
               for(190 I=1,30
               for(190 J=1,20
             190NDIST(I,J,K)=0
               for(200 I=1,30
               for(200 J=1,20
             200NDIST(I,J,6)=0
                REWIND 2
             210READ (2,END=230) INDX,MODE,JC,JF,MP,IE,IP,IZC,INC,FPROB
                IE1=IE/IB+1
                IP=(IP+9)/10
                GFS=GFS+FPROB
                SFE=SFE+IE*FPROB
                TFE=TFE+IE**2*FPROB
                SFS=SFS+JC*FPROB
                TFS=TFS+JC**2*FPROB
                SFM=SFM+FPROB*MP
                TFM=TFM+FPROB*MP**2
//              CHANGE SCALE TO 1 MEV BINS
                if(IE1>30) IE1=30
                JC1=JC/5+1
                if(JC1>20) JC1=20
                PF(IE1)=PF(IE1)+FPROB
                NPF(IE1)=NPF(IE1)+1
                if(MODE==6)THEN
                KFISS(IZC,INC)++;
                YJM(IE1,JC1)=YJM(IE1,JC1)+fabs(MP)
                goto 220
               }//endIF
//              if(MODE==6) goto 220
                if(IZCS(INDX)<0) goto 220
                NDIST(IE1,JC1,MODE)=NDIST(IE1,JC1,MODE)+1
//              TRACEBACK NDIST IS RESIDUES ONLY
                NTOT(MODE)++;
                ECM(MODE)+=IP;
                DSPIN(MODE)+=JC-JF;
                goto 210
             220 NTFISS(IE1,JC1,MODE)++;
                if(MODE==6) goto 210
                NTOTF(MODE)++;
                ECMF(MODE)+=IP;
                DSPINF(MODE)+=JC-JF
                goto 210
             230WRITE(10, 570)
                NZN=0;
               for(240 I=1,NCASC
                if(IZCS[I]>=0) NZN++;
             240}
               for(260 I=1,30
               for(260 J=1,20
                SUM=NTFISS(I,J,6)+1.E-19
                SUMM=SUM
               for(250 K=1,5
                SUM=SUM+NDIST(I,J,K)+NTFISS(I,J,K)
             250}
//              YJM(I,J)=YJM(I,J)/SUM
                YJM(I,J)=YJM(I,J)/SUMM
             260}
                if(IDIST==0) goto 280
                WRITE(10, 440 )
               for(270 I=1,30
                IEE1=IEBIN[I]-IB
                WRITE(10, 450)IEE1,IEBIN[I],(YJM(I,J),J=1,17)
             270}
                WRITE(10, 560) (JBIN[I],I=1,16)
             280WRITE(10, 460)
               for(320 K=1,5
                if(NTOT(K)==0) goto 320
                AM=NTOT(K)
                AM=AM/NZN
                EC=ECM(K)/NTOT(K)-.5
                DSP=DSPIN(K)/NTOT(K)
                WRITE(10, 580) NTITL(K),NTOT(K),NZN,AM,EC,DSP
                WRITE(10, 90001)
               for(300 I=1,30
                ISUM=0
                IAVR=0
                ISGM=0
               for(290 J=1,20
                IAVR=IAVR+NDIST(I,J,K)*(5*J-3)
                ISGM=ISGM+NDIST(I,J,K)*(5*J-3)**2
                ISUM=ISUM+NDIST(I,J,K)
             290}
                if(ISUM==0) goto 300
                AVR=IAVR
                AVR=AVR/ISUM
                STD=ISGM
                STD=Zsqrt(STD/ISUM-AVR**2)
                IEE1=IEBIN[I]-IB
                WRITE(10, 530) IEE1,IEBIN[I],(NDIST(I,J,K),J=1,16),ISUM,AVR,STD
             300}
               for(310 J=1,20
                JSUM[J]=0
               for(310 I=1,30
             310JSUM[J]=JSUM[J]+NDIST(I,J,K)
                WRITE(10, 550) JSUM
             320}
                WRITE(10, 590)
                SUMPF=0.
               for(330 I=1,30
                NI=31-I
                if(NPF(NI)<=0) goto 330
                SUMPF=SUMPF+PF(NI)
                PF(NI)=PF(NI)/NPF(NI)
                IEE1=IEBIN(NI)-IB
                WRITE(10, 600) IEE1,IEBIN(NI),PF(NI)
             330}
                SUMPF=SUMPF/NCASC
                WRITE(10, 610) SUMPF
                if(SUMPF<1.E-6) return
                SFE=SFE/GFS
                TFE=Zsqrt(TFE/GFS-SFE**2)*2.36
                SFS=SFS/GFS
                TFS=Zsqrt(TFS/GFS-SFS**2)*2.36
                SFM=SFM/GFS
                TFM=sqrt(TFM/GFS)
                WRITE(10, 640) SFE,TFE,SFS,TFS
                if(IDIST>0) WRITE(10, 650) SFM,TFM
                NTT=0
               for(340 I=1,30
               for(340 J=1,20
               for(340 K=1,6
             340NDIST(I,J,K)=NTFISS(I,J,K)
               for(350 I=1,30
               for(350 J=1,20
             350NTT=NTT+NDIST(I,J,6)
                if(NTT==0) return
                WRITE(10, 620) NTT
               for(420 K=1,6
                if(K==6) goto 360
                if(NTOTF(K)==0) goto 420
                AM=NTOTF(K)
                AM=AM/NTT
                EC=ECMF(K)/NTOTF(K)-.5
                DSP=DSPINF(K)/NTOTF(K)
                WRITE(10, 430) NTITL(K),NTOTF(K),AM,EC,DSP
                goto 380
             360WRITE(10, 630) NTITL(K),NTT
                SE=0.
                TE=0.
                SS=0.
                TS=0.
               for(370 I=1,30
               for(370 J=1,20
                SPIN=5*J-3
                EXEN=IEBIN[I]-.5*IB
                SE=SE+EXEN*NDIST(I,J,6)
                TE=TE+EXEN**2*NDIST(I,J,6)
                SS=SS+SPIN*NDIST(I,J,6)
                TS=TS+SPIN**2*NDIST(I,J,6)
             370}
                SE=SE/NTT
                TE=Zsqrt(TE/NTT-SE**2)*2.36
                SS=SS/NTT
                TS=Zsqrt(TS/NTT-SS**2)*2.36
                WRITE(10, 640) SE,TE,SS,TS
             380}
                WRITE(10, 90001)
               for(400 I=1,30
                ISUM=0
                IAVR=0
                ISGM=0
               for(390 J=1,20
                IAVR=IAVR+NDIST(I,J,K)*(5*J-3)
                ISGM=ISGM+NDIST(I,J,K)*(5*J-3)**2
             390ISUM=ISUM+NDIST(I,J,K)
                if(ISUM==0) goto 400
                AVR=IAVR
                AVR=AVR/ISUM
                STD=ISGM
                STD=STD/ISUM-AVR**2
                STD=Zsqrt(STD)
                IEE1=IEBIN[I]-IB
                WRITE(10, 530) IEE1,IEBIN[I],(NDIST(I,J,K),J=1,16),ISUM,AVR,STD
             400}
               for(410 J=1,20
                JSUM[J]=0
               for(410 I=1,30
             410JSUM[J]=JSUM[J]+NDIST(I,J,K)
                WRITE(10, 550) JSUM
             420}
//
                WRITE(10,570)
                WRITE(10,1990)
  1          990FORMAT(/' ******** DETAILS OF MULTI-CHANCE FISSION *******'
                1       '*'//,32X,'Z',6X,'N',17X,'YIELD'/)
               for(1988 I=1,110
               for(1988 J=1,180
                if(KFISS(I,J)>0)THEN
                WRITE(10,1987)I,J,KFISS(I,J)
  1          987FORMAT(30X,I3,5X,I3,10X,I8)
               }//endIF
  1          988}
*         ***** ******************************************************************
                REWIND 2
  8          102READ (2,END=8100) INDX,MODE,JC,JF,MP,IE,IP,IZC,INC,FPROB
                if(MODE==4 || MODE==5)MODE=4
                if(MODE==6)MODE=5
//          the maximum number of cascade should be less than 10
               for(8101 J=1,10
                if(ICASC(INDX,J)==0)THEN
                ICASC(INDX,J)=10**(MODE-1)
                if(ICASC(INDX,13)<JC)ICASC(INDX,13)=JC
                if(ICASC(INDX,14)>JC)ICASC(INDX,14)=JC
                if(MODE==3)THEN
                ICASC(INDX,15)=IE
               }//endIF
                goto 8102
               }//endIF
  8          101}
                goto 8102
  8          100DO 8103 I=1,NCASC
               for(8103 J=1,10
  8          103ICASC(I,11)=ICASC(I,11)+ICASC(I,J)
               for(8104 I=1,NCASC
                if(ICASC(I,11)<=0)goto 8104
               for(8104 J=I+1,NCASC
                if(ICASC(I,11)==ICASC(J,11))THEN
               for(8211 K=1,10
  8          211if(ICASC(I,K)!=ICASC(J,K))goto 8104
                ICASC(I,12)=ICASC(I,12)+1
                if(ICASC(I,13)<ICASC(J,13))ICASC(I,13)=ICASC(J,13)
                if(ICASC(I,14)>ICASC(J,13))ICASC(I,14)=ICASC(J,13)
                ICASC(J,11)=-ICASC(J,11)
               }//endIF
  8          104}
                WRITE(10,8105)
  8          105FORMAT(//' ****************** CASCADE SEQUENCE **************
                1       10X,' CASCADE',5X,'YIELD',5X,'SPIN REGION',4X,
                1       ' Eex FOR ALPHA DECAY'//)
               for(8106 I=1,NCASC
                if(ICASC(I,11)<=0)goto 8106
               for(8107 J=1,10
                if(ICASC(I,J)==1)JCASC(I,J)='N'
                if(ICASC(I,J)==10)JCASC(I,J)='P'
                if(ICASC(I,J)==100)JCASC(I,J)='A'
                if(ICASC(I,J)==1000)JCASC(I,J)='G'
  8          107if(ICASC(I,J)==10000)JCASC(I,J)='F'
                WRITE(10,8108)(JCASC(I,K),K=1,10),ICASC(I,12)+1,ICASC(I,14),
                1       ICASC(I,13),ICASC(I,15)
  8          108FORMAT(10X,10A1,3X,I6,5X,I3,' - ',I3,8X,I5)
  8          106}
                return
//
             430FORMAT (//' PRE-FISS DECAY.  MODE = ',A4,4X,'TOTAL = ',I4,4X,'MULT
               1IPLICITY =',F5.2,5X,'AVERAGE ECM =',F6.2,5X,'AVERAGE SPIN REMOVED
               2=',F5.1/1X,35(1H*)//)
             440FORMAT ('1 ABS VALUE OF M STATES AT INTERMEDIATE J VS EX FOR',
                1       ' FISSION '//)
             450FORMAT (1X,I3,'-',I3,3X,17F6.1)
             460FORMAT (' EVAPORATION RESIDUE DECAY SUMMARY'/' *******************
               1**************'//)
             470FORMAT ('1',/////' ********* COMPLETE TRACEBACK DIAGNOSTIC OF PART
               1ICLE AND GAMMA EMISSION *********'///)
             480FORMAT (/////20X,'COMPONENTS  OF  GAMMA  SPECTRUM *** '//16X,' ENE
               1RGY        E1-SPEC  E2-SPEC'//(14X,F5.2,' -',F5.2,5X,I5,5X,I5))
             490FORMAT (14X,'ABOVE  ',F5.2,5X,I5,5X,I5)
             500FORMAT (1H1)
             510FORMAT (///' Z N =',2I4,'   MODE = ',A4,'   TOTAL = ',I6,'  OUT OF
               1',I4,' EVENTS.    MULTIPLICITY = ',F5.2,'  AVERAGE ECM = ',F5.2,'
               2MEV'/1X,13(1H*),3X,'AVERAGE SPIN REMOVED = ',F5.1//)
             520FORMAT (' *************  FISSION FOLLOWS  PARTICLE ','EMISSION CAS
               1CADES  ************* '//)
             530FORMAT (1X,I3,'-',I3,3X,17I6,2F6.1)
             540FORMAT (/'  EX / J',I9,15I6,'   SUM  AVRG  STDV'/)
             550FORMAT (' SUM',7X,20I6//)
             560FORMAT (/'  EX / J',I9,15I6/)
             570FORMAT (1H1)
             580FORMAT (///' DECAY TO E.R.-S   MODE = ',A4,'   TOTAL = ',I6,'  OUT
               1OF ',I4,' EVENTS.    MULTIPLICITY = ',F5.2,'  AVERAGE ECM = ',F5.2
               2,' MEV'/1X,13(1H*),3X,'AVERAGE SPIN REMOVED = ',F5.1//)
             590FORMAT ('1 FISSION PROBABILITY AS FUNCTION OF EXCITATION'//'
               1EX. EN.    PROBABILITY     '/)
             600FORMAT (17X,I3,'-',I3,6X,1PE10.2)
             610FORMAT (' ***** TOTAL SUM OF FISSION PROBABILITIES ',1PE11.2)
             620FORMAT (///'1FISSION SUMMARY. TOTAL EVENTS = ',I4/' **************
               1**********************'///' PARICLES PRECEEDING FISSION'//' ******
               2********************* '//)
             630FORMAT (///' MODE = ',A4,4X,' TOTAL NUMBER = ',I4/1X,35(1H*)//)
             640FORMAT (//' EXCITATION ENERGY WINDOW - AVERAGE = ',F6.1,'
               1FWHM = ',F5.1//' SPIN WINDOW              - AVERAGE = ',F6.1,'
               2FWHM = ',F5.1//)
             650FORMAT (' AVERAGE  SPIN PROJECTION ',F5.1,'       AVERAGE RMS  SPI
               1N PROJECTION ',F5.1//)
               }//end
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
*/
char* LMNT (int IZ)
{
char *IELEM[103] = { "  ","H ","HE","LI","BE","B ","C ","N ","O ","F ","NE","NA"
               ,"MG","AL","SI","P ","S ","CL","AR","K ","CA","SC","TI","V ","CR"
               ,"MN","FE","CO","NI","CU","ZN","GA","GE","AS","SE","BR","KR","RB"
               ,"SR","Y ","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN"
               ,"SB","TE","I ","XE","CS","BA","LA","CE","PR","ND","PM","SM","EU"
               ,"GD","TB","DY","HO","ER","TM","YB","LU","HF","TA","W ","RE","OS"
               ,"IR","PT","AU","HG","TL","PB","BI","PO","AT","RN","FR","RA","AC"
               ,"TH","PA","U ","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO" };

IZ=max(IZ,0);
if(IZ>102) IZ=0;
return IELEM[IZ];
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*
void TCCAL (int JAH, int JZ, int JA, int NE, int KUC, int JMT,
            double BTT[6],double SIV,double BMV,double BV[5],double BW[5],double BWV[5],
            double SCX,double  EP[37], double FTHRS, double FCUT,double BSO[5],
            double SSOT, double BSOT, int IRAD, int JMESH, double TP[37][81], int LMX[37])
{
//              double TCCAL PACKAGE , E.D. ARTHUR, T-2, L.A.N.L. CONTAINS
//              FOLLOWING -
//                       PRPARM
//              ZEROS ->(  SERIES    INT1      INT2    INT3    INT4  )
//              JMTRAN ->( RPOT  NEUTS   COOL   OUTC    FAZE2
//
//              ADAPTED FOR VAX BY A. GAVRON, 3/28/84
//
//              PACKAGE CALCULATES TRANSMISSION COEFFICIENTS
//
//              CALCULATES PROBABILITIES PARTICLE   EMISSION CHANNELS
//

//double P[81], ATT[6], AV[5], AW[5] ;

//                COMMON /Q55/ ASO(4)
//                COMMON /Q1/ EN, EC, IZ, IA, IS, IMDL, SI, SII, A, FACT, FACT1,IPUN
//                COMMON /Q3/ IAH
//                COMMON /Q30/ AW1(4), W1, ETHRS, ECUT
//                COMMON /Q99/ RSOT, ASOT, RIV, AIV, AWV(4), WVO

char *XPART1[7] = {"   ","NEUT","PROT","ALPH","DEUT","TRIT","HE-3"};

double  PC1[7]={0.,0.,1.,2.,1.,1.,2.};
double  PM1[7]={0.,1.00866522,1.007825,4.002635,2.0141,3.01605,3.016035};

//
                IZ=JZ;
                IMT=JMT;
                RIV=SIV;
                int IAH=JAH;
                JUC=KUC;
                ETHRS=FTHRS;
                IMESH=JMESH;
                RSOT=SSOT;
                int IA=JA;
                RCX=SCX;
                ECUT=FCUT;
                AIV=BMV;
                ASOT=BSOT;
               for(int I=1; I<=4; I++) {
                  AV[I]=BV[I];
                  AW[I]=BW[I];
                  ASO[I]=BSO[I];
                  AWV[I]=BWV[I];
                  }

               for(int I=1; I<=5; I++)  ATT[I]=BTT[I];
//
//              PARTICLE    IAH
//              ********    ***
//              NEUTRONS    1
//              PROTONS     2
//              ALPHAS      3
//              DEUTERONS   4
//              TRITONS     5
//              HE3         6
//
//              CALCULATE CONTINUUM CHANNEL PROBABILITIES
                TC=IZ;
                IIM=IMT;
                TM=IA;
                fprintf(f04,"%s particle emission\n",XPART1[IAH]);
                PC=PC1[IAH];
                PM=PM1[IAH];
                int NEP=0;
                int JUP2=2*JUC;
//
//              ECUT EQ 0 IGNORE
//              FIRST SET OF IMAG PARAM FOR IMAG SURFACE TERM
                if(NE==0) goto 70;
//
//              SET UP CONTINUUM CALCULATION
                NOP=JUP2;
//
//              STORE ZERO TRANSMISSION COEFFICIENTS FOR ZERO INCIDENT ENERGY
//
//              SET UP OPTICAL MODEL PARAMETERS
                TC=(1.0-PCS)*TC;
               for(int I=1; I<=4; I++) AT[I]=ATT[I];

                PRPARM (AV,AW);
//
                WVOLF=0.
//              STORE TRANSMISSION COEFFICIENTS
               for(60 N=1,NE
//
//              CALCULATE COEFFICIENTS
                V0=AV[1]+AV[2]*EP[N]+AV[3]*EP[N]*EP[N]+AV[4]*log(EP[N])
                W0=AW[1]+AW[2]*EP[N]+AW[3]*EP[N]*EP[N]+AW[4]*log(EP[N])
                ATT[5]=ASO[1]+ASO[2]*EP[N]+ASO[3]*EP[N]*EP[N]+ASO[4]*log(EP[N])
                AT[5]=ATT[5]
                WVO=AWV[1]+AWV[2]*EP[N]+AWV[3]*(EP[N]*EP[N])+AWV[4]*log(EP[N])
                WVO=max(WVO,0);
                V0=max(V0,0);
                W0=max(W0,0);
                if(ECUT>0.) W1=AW1[1]+AW1[2]*(EP[N]-ETHRS)+AW1[3]*(EP[N]*EP[N])+AW1[4]*log(EP[N])
                if(ECUT>0.&&EP[N]>ECUT) W0=max(W1,0.)
                VS=AT[5]
                VS=max(VS,0.)
                JMTRAN (P,EP[N],LMAX,IRAD)
               for(40 I=1,NOP
                TP(N,I)=0.
                if(I>LMAX) goto 40
                TP(N,I)=P[I]
              40}
                AMU=PM*TM/(PM+TM)
                ECM=EP[N]
                ALAMDA=4.651/sqrt(AMU*ECM)
                SUM=0.
               for(int I=1; I<=LMAX; I++) {
                        AL=I-1;
                        SUM+=(2.*AL+1.)*P[I];
                        }

                SUM*=31.4159*ALAMDA*ALAMDA;
                LMX[N]=LMAX;
//              WRITE(4,70) EP[N],ALAMDA,SUM,LMX[N]
//           TC WRITE (4,80) (TP(N,I),I=1,LMAX)
              60}
//
//              PUNCH COEFFICINTS IF REQUESTED
//
              70}
                return
//
//         TC80 FORMAT (6E16.5)
//
               }

*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
/*                double PRPARM (AV,AW)
//
//              PRINTS SPHERICAL OPTICAL MODEL PARAMETERS
//
//              AV - ARITHMETIC EXPRESSION FOR THE REAL WELL DEPTH
//              AW - ARITHMETIC EXPRESSION FOR THE IMAGINARY WELL DEPTH
//
                DIMENSION AV(4), AW(4), SOG(2,3)
                COMMON /Q30/ AW1(4), W1, ETHRS, ECUT
                COMMON /Q55/ ASO(4)
                COMMON /Q99/ RSOT, ASOT, RIV, AIV, AWV(4), WVO
//              DATA SOG /60HSAXON DERIVATIVE            GAUSSIAN        VOLUME SA
//             1XON    /
//
                WRITE (4,50)
                WRITE (4,60)
                WRITE (4,70) PC,PM,TC,TM
                if(IIM!=4)goto 20
              10WRITE (4,80) AV,AW,AW1
                goto 30
              20}
//              WRITE (4,70) (SOG(I,IIM),I=1,2)
                WRITE (4,40)
                WRITE (4,90) AV
                WRITE (4,100) AW
                if(ECUT>0.) WRITE (4,110) ECUT,AW1
              30}
                WRITE (4,120) (AT(K),K=1,4),ASO
                WRITE (4,130) RSOT,ASOT
                WRITE (4,140) (AWV(K),K=1,4),RIV,AIV
//
                return
//
              40FORMAT (20X,"SAXON REAL WELL",26X,"IMAGINARY WELL"/)
              50FORMAT (/)
              60FORMAT (10X,"TRANSMISSION COEFFICIENTS CALCULATED FROM THE ","FOLL
               1OWING"/10X,"OPTICAL MODEL PARAMETERS"/)
              70FORMAT (20X,'PROJECTILE     CHARGE',F5.1,10X,'MASS',1PE14.6/20X,'T
               1ARGET       CHARGE',0PF5.1,10X,'MASS',1PE14.6/)
              80FORMAT (20X' SAXON REAL WELL AND COMBINATION VOL AND SAXON DERIVAT
               1IVE IMAG WELL'/20X' V = '4F8.2/20X' W(VOL) = '4F8.2/20X' W(D) = '
               24F8.2)
              90FORMAT (20X' V = 'F8.3' + ('F8.3') X E + ('F8.3') X E2 + ('F8.3')
               1X LN(E)'/)
             100FORMAT (20X' W = 'F8.3' + ('F8.3') X E + ('F8.3') X E2 + ('F8.3')
               1X LN(E)'/)
             110FORMAT (20X'ABOVE E =  'F10.3'  MEV,THE FOLLOWING SURF DERIV POT I
               1S USED'/20X' W = 'F8.3' + ('F8.3' )XE + ('F8.3' )XE2 +( 'F8.3' )XL
               2N(E)'/)
             120FORMAT (20X,'RR = ',1PE14.6/20X,'AR = ',1PE14.6/20X,'RI = ',1PE14.
               16/20X,'AI = ',1PE14.6//20X,'VSO = ',0P,F8.3' + ('F8.3') X E + ('F
               28.3')X E2 + ('F8.3') X LN(E)'/)
             130FORMAT (20X'RSO = ',1PE14.6/20X'ASO = '1PE14.6//)
             140FORMAT (/20X' THE FOLLOWING SAXON VOL PARAM ARE USED '/20X'WVO= 'F
               18.3' + ('F8.3')X E +('F8.3')X E2 + ('F8.3') X LN(E)'/20X'RIV = '1
               2PE14.6/20X'AIV = '1PE14.6//)
               }//end*/
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double GETSEED(void)  //  return ISEED
{        time_t t;
        int  Y=time(NULL);
        int  ISEED=Y*2+1;
        return ISEED;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double CLOCK() //return I
{                time_t t;
                int  I=time(NULL);
                return I;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double Zsqrt (double X) {return sqrt(max(X,0));  }
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int ICOUPL (double SC, double SE) { return (SE+SC==0. ? -2 : -1 );}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

