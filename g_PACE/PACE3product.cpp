#include <QTextStream>
#include <QString>
#include "WinLise_Constant.h"
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "ftype.h"
#include "ParticleStates.h"

#include <cmath>

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))

void PRODCT(int II, int IZ, int IN, int NCASC, fusion_event **a, const char *filename_cs, QTextStream& s);
void OUTRES(int ICTRL, int NP, fusion_event *f, QTextStream &s);
int  ISUM(int I, int J, int L1, int L2) ;
//extern char* ElementName(int IZ) ;
extern QString ElementName(int IZ);
extern char *s100m;
extern char *s40star;
//extern char *s35x;
const char s35x[10] = "         ";
extern FILE *f09;
extern double mzsqrt(double X);
const char *s10m="           ";

//common /oleg_N/
int NR[17][33][38];
int NRW[17][5][38];       // from SRCH

QString qbuf;
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PRODCT(int II, int IZ, int IN, int NCASC, fusion_event **a, const char *filename_cs, QTextStream& s)
{
//
//     ORIGINAL JULIAN SUBROUTINE
//
//     ***********
//
//      SORTING ROUTINE OF RESIDUAL NUCLEI
//

//=========      common/OUT/
//extern double _ERGCS[N_fusion_event];
//extern int    _IZCS[N_fusion_event],_INCS[N_fusion_event],_JCS[N_fusion_event],_MJCS[N_fusion_event];
//extern double _VY[N_fusion_event],_VZ[N_fusion_event],_VX[N_fusion_event];

//=========      common /XQSIG/
extern double _SIGMA, _ENERGY, _LowLimit, _HighLimit;
//extern int _NCASC;
//extern double _DSIG[19];
extern double _AJNUC;
extern int _INPUT;

//=========      common/XQANG/
//extern double _SUM[11][301],_CT;     // SUM - transposed
extern int _MDIR;
extern int _IDIST;

//=============================================================
const char *MA1;
double PRES[8006];
int    NRES[8006], /*NNUC[8006],*/ IZNUC[8006] ,INNUC[8006];
bool flagA6=false;


int J, IRES, IFISS, NC,I, IZC, INC, IA, IBIN, IAT, K,NRESS;
double SRES,SIGS,PRESS, SIG, PR;
int  IBIN1,IBIN2,IZR, INR, /* ID, */ JK, IK;

//=============================================================

FILE *file_cs;
file_cs=fopen(filename_cs,"wt");


if(file_cs) {
   fprintf(file_cs,"!-title- Zcomp  Ncomp  Energy  CSfus  Input _AJNUC\n");
   fprintf(file_cs,"!PACE4 %d %d %.3f %.3g %d %.3g",IZ,IN,_ENERGY,_SIGMA, _INPUT, (_INPUT==1 ? 0: _AJNUC));
}

for(J=1; J<=8005; J++) NRES[J]=0;

IRES=0;
IFISS=0;
NC=NCASC+II-1;

//--------------------------------------------------------
for(I=II; I<=NC; I++)
    if(a[I]->Z<0) IFISS++;
    else  {
        IZC=IZ-a[I]->Z;
        INC=IN-a[I]->N;
        IA=IZC+INC;
        IBIN=IA*(IA+1)/2+IZC+1;

        if(IBIN>8005) {
            if(IRES!=1) {
                fprintf(f09,"\\par THE FOLLOWING PRODUCTS WERE NOT BINNED\\par   Z    N");
                s << "<p> THE FOLLOWING PRODUCTS WERE NOT BINNED</p><p>   Z   N  </p>";
                IRES=1;
            }
            fprintf(f09,"\\par %3d %4d",a[I]->Z,a[I]->N);
            s << "<p> " << qbuf.number(a[I]->Z) << " " << qbuf.number(a[I]->N) << "</p>";
        }
        else    NRES[IBIN]++;
    }
//--------------------------------------------------------

fprintf(f09,"\\par\\par");
s << "<p>&nbsp;</p>";
fprintf(f09,"\\par\\b\\cf1\\fs24        1.Yields of residual nuclei\\fs16\\cf0\\b0"
            "\\par\\par\\b\\cf4\\fs20     Z   N    A        events     percent    x-section(mb)\\b0\\cf0");
s <<"<h2 align=center>1. Yields of residual nuclei</h2><table align=center cellpadding=10><tr style=\"color:green\"><th>Z</th><th>N</th><th colspan=2>A</th><th>events</th><th>percent</th><th>x-section(mb)</th></tr>";

      SRES=NCASC;
      IAT=IZ+IN;
      K=1;
      SIGS=0.;
      NRESS=0;
      PRESS=0.;
//--------------------------------------------------------
    for(I=1; I<=89; I++) {
                        //     IF DIMENSION CHANGED FROM 8005, CHANGE THIS      DO RANGE.
              IBIN1=I*(I-1)/2+1;
              IBIN2=I*(I+1)/2;
              IZR=IZ;
              IAT=IZ+IN-I+1;

              for(J=IBIN1; J<=IBIN2; J++) {

                  if(NRES[J]!=0) {
                      PRES[J]=double(NRES[J])/SRES*100.;
                      INR=IAT-IZR;
                      SIG=_SIGMA*PRES[J]/100.;

                      double SIG_ER = SIG / sqrt(NRES[J]);

                      MA1=ElementName(IZR).toStdString().c_str();
                      QString testName = ElementName(IZR);
                     // qDebug()  << testName;
                      fprintf(f09,"\\par    %3d %3d\\b   %3d %.2s\\b0     %6d    %7.3g%%  %9.3g",
                                        IZR,INR,IAT,MA1,NRES[J],PRES[J],SIG);
                      s << "<tr><td align=center>" << qbuf.number(IZR) << "</td><td align=center>" << qbuf.number(INR) << "</td><td><b>" << qbuf.number(IAT) << "</b></td><td><b>" << testName << "</b></td><td align=center>"  << qbuf.number(NRES[J])
                           <<  "</td><td align=center>" << qbuf.number(PRES[J],'g',3) << "&#37;</td><td align=center>" << qbuf.number(SIG,'g',3) << "</td></tr>";
                      if(IAT<6) flagA6=true;

                      if(file_cs) fprintf(file_cs,"\n%d %d %10.3g %10.3g",IZR,INR,SIG,SIG_ER);


                      if(PRES[J]>=_LowLimit && K<16 && _HighLimit >=PRES[J] ) {
                        IZNUC[K]=IZR;
                        INNUC[K]=INR;
                       //Qt-Oleg NNUC[K]=NRES[J];
                        K++;
                        }

                      SIGS+=SIG;
                      PRESS+=PRES[J];
                      NRESS+=NRES[J];
                      }
                   IZR--;
                   }
              }

//--------------------------------------------------------

      PR=double(IFISS)/SRES*100.;
      SIG=_SIGMA*PR/100.;
      PRESS+=PR;
      NRESS+=IFISS;
      if(IFISS>0){
          fprintf(f09,"\\par\\cf3\\b     Total fission    \\b0  %6d    %7.3g%%  %9.3g\\cf0",IFISS,PR,SIG);
          s << "<tr><td colspan=3><b>  Total fission</b>   </td><td></td><td align=center>" << qbuf.number(IFISS) << "</td><td align=center>" << qbuf.number(PR) << "</td><td align=center> " << qbuf.number(SIG) << "</td></tr>";
      }
      SIGS+=SIG;
      fprintf(f09,"\\par\\cf4     \\b TOTAL             %6d    %7.3g%%  %9.3g\\b0\\fs16\\cf0",NRESS,PRESS,SIGS);
      s << "<tr style=\"color:green\"><td colspan=2>  <b>TOTAL </b> </td><td></td><td></td><td align=center>" << qbuf.number(NRESS) << "</td><td align=center>" << qbuf.number(PRESS) << "</td><td align=center> " << qbuf.number(SIGS) << "</td></tr>";
      fprintf(f09,"\\par");
      s << "</table>";
      fclose(file_cs);

      if(flagA6) {
          fprintf(f09,"\\par\\par\\fs12     The Cascade is stopped if a mass of residual less than 6 \\fs16 ");
          s << "<p>    The Cascade is stopped if a mass of residual less than 6 </p>";
      }

      if(_IDIST==0)return;

      OUTRES(1,0,a[1],s);
      K--;

//========================================================
  for(I=II; I<=NC; I++)
      {

      if(a[I]->Z<0)  continue;

      for(J=1; J<=K; J++)
         if(a[I]->Z == IZNUC[J] &&
            a[I]->N == INNUC[J]) break;

      if(J>K)   JK=16;
      else     { a[I]->MJ=J;  JK=J; }

      OUTRES(2,JK,a[I],s);
      }
//========================================================

      K++;

      fprintf(f09,"\\par\\par\\par\\b\\cf1\\fs24       2.Angular distribution results\\fs16\\cf0\\b0");
      s << "<p>&nbsp;</p><h2 align=center>  2. Angular distribution results</h2>" ;
      if(K>16)K=16;
      fprintf(f09,"\\par\\par\\cf2\\f2\\fs20  ");
      if(_MDIR==0){
          fprintf(f09,"*** Spin alignment perpendicular to recoil axis  - standard compound nucleus angular distribution");
          s << "<p>*** Spin alignment perpendicular to recoil axis  - standard compound nucleus angular distribution</p>";
      }
      if(_MDIR==1){
          fprintf(f09,"*** Spin alignment perpendicular to reaction plane  - angular distribution is around Z axis perpendicular to"
                                "\\par  *** reaction plane");
          s << "<p>*** Spin alignment perpendicular to reaction plane  - angular distribution is around Z axis perpendicular to <br> *** reaction plane</p>";
      }
      fprintf(f09,"\\cf0\\f0\\fs16  ");

  for(I=1; I<=K; I++) {
      fprintf(f09,"\\cf3\\f2\\fs20\\b \\par\\par%s     2.%d  Energy and angular distribution of ",s35x,I);
      s << "<br><h3>" <<  "2." << qbuf.number(I) << " Energy and angular distribution of ";
      if(I==K) {fprintf(f09,"\\cf1  ALL \\cf3  residual nuclei"); s << "ALL residual nuclei</h3>";}
      else     {
          fprintf(f09,"residual nucleus  Z=\\cf1 %d\\cf3    N=\\cf1 %d\\cf3   \\cf2 (%d%.2s)\\cf3 ",
                                IZNUC[I],INNUC[I],IZNUC[I]+INNUC[I],ElementName(IZNUC[I]).toStdString().c_str());
          s << " residual nucleus Z = <span style=\"color:blue\">" << qbuf.number(IZNUC[I]) << "</span> and  N = <span style=\"color:blue\">" << qbuf.number(INNUC[I]) << "</span> ("
               << qbuf.number(IZNUC[I]+INNUC[I]) << ElementName(IZNUC[I])[0] <<  ElementName(IZNUC[I])[1] << ")</h3>";
      }
      fprintf(f09,"\\cf0\\f0\\fs14\\b0 ");

//      fprintf(f09,"\\par%s",s100m);
      IK=I;
      if(IK==K) IK=16;
      s.flush();
      OUTRES(3,IK,a[I],s);
      fprintf(f09,"\\fs16 ");
  }
  s.flush();
  return;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void OUTRES(int ICTRL, int NP, fusion_event *f, QTextStream& s)
{

//
//     ORIGINAL JULIAN SUBROUTINE
//
//     ***********
//
//   ICTRL IS CONTROL PARAMETER
//   ERESL IS KINETIC ENERGY OF RESIDUAL NUCLEUS IN LABORATORY
//   ARESL IS LABORATORY ANGLE OF RESIDUAL NUCLEUS
//

extern int _INPUT;
//=========      common /XQSIG/
extern double _SIGMA;
extern int _NCASC;
//extern double _DSIG[19];
//=========      common/XQANG/
//extern double _SUM[11][301],_CT;     // SUM - transposed
extern int    _MDIR;
//=========      common/OUT/

//==========================
//=========      common/SRCH/
//extern double _PROB[5556],_EMAX[5],_EMIN[5],_BE[5],_GE1,_GE2;
//extern int    _IPROB,_MJ[5556],_IMODL;


//=========      common/XQRES/
extern double _ENERGY,_VZC;
extern int    _IZ,_IN;//,_ITRAC;


//==========================
double DSIG[37];
double EREC;
int ILOW, IDELE, IREC, K, I, L, M, KW,IS,IJ,LQ;

double AM;
double A, E2, AK, E1, FAC, TET;

static int IWR=0;
static double ELOW, DELE,DELANG;
//static double ANGLE1[37],ANGLE2[37];
static double ANGLE1[38],ANGLE2[38];//MPK fixing out of bound arrays
//static double EWR[4],EWR1[4],EWR2[4];
static double EWR[5],EWR1[5],EWR2[5]; //MPK fixing out of bound arrays
//
//   NR(ENERGY,ANGLE) IS DISTRIBUTION OF RESIDUAL NUCLEI
//   NRW(BIN,ANGLE) IS DISTRIBUTION OF RESIDUAL NUCLEI IN 4 ENERGY BINS
//   EWR  ARE SEPARATION ENERGIES BETWEEN FOUR ENERGY BINS
//   EWR1 LOW ENERGIES OF FOUR ENERGY BINS
//   EWR2 HIGH ENERGIES OF FOUR ENERGY BINS. FOURTH IS INFINITY
//

if(ICTRL<2)
 {                              //   INITIALIZATION
      A = _IZ+_IN;
      EREC=460.*A*_VZC*_VZC;
//      ILOW=(1.-0.4*_ENERGY/A)*EREC +1.001;
      ILOW=(1.-0.4*_ENERGY/A)*EREC +0.001;
      ELOW=max(ILOW,0);
      IDELE=(EREC-ELOW)/16.+0.5;
      DELE=max(IDELE,1);
      DELANG=1.+4.*_MDIR;
      IWR=1;
      EWR[1]=ELOW;
      IREC=EREC;
      EWR[2]=IREC;
      if(EWR[2]<EWR[1]+1.)EWR[2]=EWR[1]+1.;
      EWR[3]=2.*EWR[2]-ELOW;
      if(DELANG==0.||DELANG>5.)DELANG=5.;
//
//   ENERGY INTERVAL OF 1ST OF 30 BINS IS FROM 0. TO ELOW
//   DELE IS ENERGY INTERVAL IN OTHER 29 BINS
//      31ST BIN FOR OVERFLOW
//      32ND BIN IS FOR TOTAL
//   DELANG IS ANGLE INTERVAL IN 18 ANGLE BINS
//   IWR IS CONTROL - 0 IS NO FOUR BINS, 1 IS YES
//
//L4:   format(3F5.0,I5,3F5.0);

   for(K=1; K<=32; K++)
        for(I=1; I<=16; I++)
           for(L=1; L<=37; L++)
                        NR[I][K][L]=0;


   for(M=1; M<=36; M++) {
      AM=M;
      ANGLE1[M]=(AM-1.0)*DELANG;
      ANGLE2[M]=AM*DELANG;
      }

   if(IWR==0)return;

   for(I=1; I<=16; I++)
     for(KW=1; KW<=4; KW++)
       for(L=1; L<=37; L++)
              NRW[I][KW][L]=0;

      EWR1[1]=0.;
      EWR2[1]=EWR[1];
      EWR1[2]=EWR[1];
      EWR2[2]=EWR[2];
      EWR1[3]=EWR[2];
      EWR2[3]=EWR[3];
      EWR1[4]=EWR[3];
}
//==============================================
//
//   MONITORING
//
else if(ICTRL==2) {

      if(f->y<ELOW) K=1;
      else              K=min(31, (f->y-ELOW)/DELE + 2);

      L = min(37,f->z/DELANG+1);

      if(NP<=15) {
              NR[NP][K][L]++;
              NR[NP][32][L]++;
              }

      NR[16][K][L]++;
      NR[16][32][L]++;

      if(IWR==0)return;


             if(f->y<=EWR[1])  KW=1;
      else   if(f->y<=EWR[2])  KW=2;
      else   if(f->y<=EWR[3])  KW=3;
      else                    KW=4;

      if(NP<=15)  NRW[NP][KW][L]++;

      NRW[16][KW][L]++;
      }
//====================================================================
//   GRAPHING
//

else if(ICTRL>2) {

int k;

char s500[134]="\\par\\fs14 "
"---------------------------------------------------------------------------------------------------------------------------";
char s501[135]="\\par "
"--------------|-----------------------------------------------------------------------------------------------------------|-----|";
char s502[129]="\\par "
"--------------|------------------------------------------------------------------------------------------------------------";
char s503[129] ="\\par "
"--------------|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|";
char s30[129] = "\\par "
" Energy range |_________________________________________ Angular_range_(deg)______________________________________________|";

//char *s1 = "\\par"
//"              |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |";

char sVelocity[106] = "\\par\\i\\fs14\\cf4         Residual Velocity/c \\i0  Vz=%.2e(sig=%.2e)    root-mean-square Vxy=%.2e\\fs16\\cf0 ";
char s32[21] = "\\par               |";
char s36[20] = "\\par %5.1f -%6.1f |";
char s38[20] = "\\par  Above %6.1f |";
char s53[21] = "\\par dSigma/dOmega |";
char s15[21] = "\\par               |";
char sT[21] =  "\\par  Total        |";
char sB[20] =  "\\par  Below %6.1f |";
char sM[21] =  "\\par    (MeV)      |";
char sD[5]="%5d|";
char sDe[7]="     |";
//char *sF="%5.0f|";
char sF1[8]="%4.0f |";
char sGb[7]="%5.1g|";
char sG[7]="   0 |";

double GlobSum=ISUM(NP, 32, 1, 37);

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//------------------------------------------------- VELOCITY if MODE=1
if(_INPUT==1 && GlobSum > 2 && NP!=16) {
  double v, vz=0, vxy, dvz, dvxy, s_vz, s_vxy, ener, angle,
          angsin, angcos,vzm, fff;

  v  = max(f->Ad(),1); // avoid fission events
  fff=2.0/(aem_MeV*v);
  s_vz = s_vxy = 0;

  //-----------------------------------------------
  for(k=1;k<=37;k++) //angle
    {
    if(k<37) angle=(ANGLE2[k]+ANGLE1[k])/2.;
    else     angle= ANGLE2[k];

    angle = (angle*RAD);
    angcos = cos(angle);

    //--------------------------------
    for(K=1; K<=31; K++)       //energy
        {
        if(NR[NP][K][k]<=0) continue;

             if(K==1) ener=ELOW;
        else if(K<31) ener=(double(K)-1.5)*DELE+ELOW;
        else          ener=30.*DELE+ELOW;

        v    = sqrt(fff*ener);
        vz   = v *angcos;

        s_vz += vz *NR[NP][K][k];
        }
    }

  vzm  = s_vz  / GlobSum;
//================================================================== deviation


  s_vz = s_vxy = 0;
  //-----------------------------------------------
  for(k=1;k<=37;k++) //angle
    {
    if(k<37) angle=(ANGLE2[k]+ANGLE1[k])/2.;
    else     angle= ANGLE2[k];

    angle = (angle*RAD);
    angsin = sin(angle);
    angcos = cos(angle);

    //--------------------------------
    for(K=1; K<=31; K++)       //energy
        {
        if(NR[NP][K][k]<=0) continue;

             if(K==1) ener=ELOW;
        else if(K<31) ener=(double(K)-1.5)*DELE+ELOW;
        else          ener=30.*DELE+ELOW;

        v    = sqrt(fff*ener);
        vz   = v *angcos;
        vxy  = v *angsin;

        s_vz  += (vz  - vzm )*(vz  - vzm )*NR[NP][K][k];
        s_vxy += (vxy)*(vxy)*NR[NP][K][k];
        }
  }

  dvz  = sqrt( s_vz  / GlobSum );
  dvxy = sqrt( s_vxy / GlobSum );

  fprintf(f09,sVelocity,vz,dvz,dvxy);
  s << "<p><em> Residual velocity/c </em> Vz = " << qbuf.number(vz,'e',2) << "(sig = " << qbuf.number(dvz,'e',2) << ") rms Vxy = " << qbuf.number(dvxy,'e',2) << "</p><p>&nbsp;</p>";
}

s.flush();

fprintf(f09,s500);
fprintf(f09,s30);
fprintf(f09,sM);
s << "<table border=1 cellspacing=0><tr><th>Energy Range</th><th colspan=18>Angular range (deg)</th></tr><tr><th>(MeV)</th>";

E2=0;
for(k=1;k<=18;k++) {
    fprintf(f09,sF1,ANGLE1[k]);
    s << "<td>" << qbuf.number(ANGLE1[k]) << "</td>";
}
s << "</tr><tr><td></td>";

//      fprintf(f09,s1);
fprintf(f09,s15);
for(k=1;k<=18;k++) {
    fprintf(f09,sF1,ANGLE2[k]);
          s << "<td>" << qbuf.number(ANGLE2[k]) << "</td>";
      }
      s << "</tr>";
      fprintf(f09,s502);

      if(ELOW>=0.01) {
          E2=ELOW;
          fprintf(f09,sB,E2);
          s << "<tr><td> Below " << qbuf.number(E2) << " </td>";
          for(k=1;k<=18;k++){
              if(NR[NP][1][k]>0) {
                  fprintf(f09,sD,NR[NP][1][k]);
                  s << "<td align=center>" << qbuf.number(NR[NP][1][k]) << "</td>";
              }else {
                  fprintf(f09,sDe);
                  s << "<td></td>";
              }
          }
          s << "</tr>";
      }


      for(K=2; K<=30; K++) {
          AK=K;
          E1=(AK-2.)*DELE+ELOW;
          E2=E1+DELE;
          if(ISUM(NP,K,1,18)>0) {
              fprintf(f09,s36,E1,E2);
              s << "<tr><td>" << qbuf.number(E1) << " - " << qbuf.number(E2) << "</td>";
              for(k=1;k<=18;k++){
                  if(NR[NP][K][k]) {
                      fprintf(f09,sD,NR[NP][K][k]);
                      s << "<td align=center>" << qbuf.number(NR[NP][K][k]) << "</td>";
                  }else  {
                      s << "<td></td>";
                      fprintf(f09,sDe);
                  }
              }
              s << "</tr>";
          }
      }

      s << "<tr><td> Above " << qbuf.number(E2) << " </td>";
      fprintf(f09,s38,E2);
      for(k=1;k<=18;k++){
                if(NR[NP][31][k]) {fprintf(f09,sD,NR[NP][31][k]); s << "<td align=center>" << qbuf.number(NR[NP][31][k]) << "</td>";}
                else     {         fprintf(f09,sDe);s << "<td></td>";}
      }
      fprintf(f09,s502);
      s << "</tr>";

      FAC=_SIGMA/(6.2832*_NCASC);
      for(I=1; I<=36; I++) {
          TET=(I-0.5)*0.01745;
          DSIG[I]=FAC*NR[NP][32][I]/(sin(TET)*0.01745);
      }

//-------------------------------- TOTAL
      fprintf(f09,sT);
      s << "<tr><td> Total </td>";
      for(k=1;k<=18;k++){
                if(NR[NP][32][k]) {
                    s << "<td align=center>" << qbuf.number(NR[NP][32][k]) << "</td>";
                    fprintf(f09,sD,NR[NP][32][k]);
                }else {
                    fprintf(f09,sDe);
                    s << "<td></td>";
                }
      }
      fprintf(f09,s53);
      s << "</tr><tr><td>dSig/dOmeg</td>";

      for(k=1;k<=18;k++){
           if(DSIG[k]>0) {fprintf(f09,sGb,DSIG[k]); s << "<td>" << qSetRealNumberPrecision(2) << (DSIG[k]) << "</td>";}
           else         { fprintf(f09,sG,DSIG[k]); s << "<td> 0.00 </td>";}
      }
      s << "</tr>";
      fprintf(f09,s503);

      if(IWR==0){return;}

      for(I=1; I<=4; I++) {
          if(I!=4){fprintf(f09,s36,EWR1[I],EWR2[I]); s << "<tr><td>" << qbuf.number(EWR1[I]) << " - " << qbuf.number(EWR2[I]) << "</td>";
          }else    {fprintf(f09,s38,EWR1[4]); s << "<tr><td> Above " << qbuf.number(EWR1[4]) << "</td>";}

          for(k=1;k<=18;k++)
              if(NRW[NP][I][k]>0){ fprintf(f09,sD,NRW[NP][I][k]); s << "<td align=center>" << qbuf.number(NRW[NP][I][k]) << "</td>";
              }else  {            fprintf(f09,sDe); s << "<td></td>";}
          s << "</tr>";
      }
      s << "</table><p>&nbsp;</p>";
      fprintf(f09,s500);
      IS=0;

      for(IJ=1; IJ<=32; IJ++)
          for(LQ=19; LQ<=37; LQ++)
              IS+=NR[NP][IJ][LQ];

      if(IS == 0) return;
      /// 18-36  degrees
      s << "<table border=1 cellspacing=0><tr><th>Energy Range</th><th colspan=18>Angular range (deg)</th></tr><tr><th>(MeV)</th>";
      fprintf(f09,s30);
      fprintf(f09,sM);

      for(k=19;k<=36;k++) {fprintf(f09,sF1,ANGLE1[k]);
          s << "<td align=center>" << qbuf.number(ANGLE1[k]) << "</td>";}
      fprintf(f09,"Above|");
      s << "</tr><tr><td> Above </td>";
      fprintf(f09,s32);
      for(k=19;k<=36;k++) {fprintf(f09,sF1,ANGLE2[k]);  s << "<td align=center>" << qbuf.number(ANGLE2[k]) << "</td>";}
      fprintf(f09,sF1,ANGLE2[36]);
      s << "</tr>";
      fprintf(f09,s501);
      if(ELOW>=0.01){
              E2=ELOW;
              fprintf(f09,sB,E2);
              s << "<tr><td>Below " << E2 <<"</td>";
              for(k=19;k<=36;k++)
                  if(NR[NP][1][k] >0) {fprintf(f09,sD,NR[NP][1][k]);s << "<td align=center>" << NR[NP][1][k] << "</td>";}
                  else             { s << "<td></td>";  fprintf(f09,sDe);}
              s << "</tr>";
      }

      for(K=2; K<=30; K++) {
      AK=K;
      E1=(AK-2.)*DELE+ELOW;
      E2=E1+DELE;
      if(ISUM(NP,K,19,37)>0) {
           fprintf(f09,s36,E1,E2);
           s << "<tr><td>" << E1 << " - " << E2 << "</td>";
           for(k=19;k<=36;k++)
                if(NR[NP][K][k]>0) {fprintf(f09,sD,NR[NP][K][k]); s << "<td align=center>" << NR[NP][K][k] << "</td>";}
                else               { fprintf(f09,sDe); s << "<td></td>";}
           s << "</tr>";
           }
      }

      fprintf(f09,s38,E2);
      s << "<tr><td> Above " << E2 << "</td>";
      for(k=19;k<=36;k++)
                if(NR[NP][31][k]>0) {fprintf(f09,sD,NR[NP][31][k]); s << "<td align=center>" << NR[NP][31][k] << "</td>";}
                else            {    fprintf(f09,sDe); s << "<td></td>";}
      s << "</tr><tr><th>Total</th>";
      fprintf(f09,s501);
      fprintf(f09,sT);
      for(k=19;k<=36;k++)
                if(NR[NP][32][k]>0) {fprintf(f09,sD,NR[NP][32][k]);s << "<th>" << NR[NP][32][k] << "</th>";}
                else               { fprintf(f09,sDe); s << "<th></th>";}

      s << "</tr><tr><td>dSig/dOmeg</td>";
      fprintf(f09,s53);
      for(k=19;k<=36;k++) 
           if(DSIG[k]>0) {fprintf(f09,sGb,DSIG[k]); s << "<td align=center>" << DSIG[k] << "</td>";}
           else         { fprintf(f09,sG,DSIG[k]); s << "<td></td>";}

      s << "</tr>";
      fprintf(f09,s501);
      if(IWR==0){s << "</table>"; return;}

     for(I=1; I<=3; I++){
         s << "<tr><td>" << EWR1[I] << " - " << EWR2[I] << "</td>";
      fprintf(f09,s36,EWR1[I],EWR2[I]);
      for(k=19;k<=36;k++)
        if(NRW[NP][I][k]>0) {fprintf(f09,sD,NRW[NP][I][k]);s << "<td align=center>" << NRW[NP][I][k] << "</td>";}
        else            { s << "<td></td>";   fprintf(f09,sDe);}
      }

      fprintf(f09,s38,EWR1[4]);
      s << "<tr><td> Above " << EWR1[4] << "</td>";
      for(k=19;k<=36;k++)
          if(NRW[NP][4][k]>0){ fprintf(f09,sD,NRW[NP][4][k]); s << "<td>" << NRW[NP][4][k] << "</td>";}
          else             {   fprintf(f09,sDe); s << "<td></td>";}

      s << "</tr></table>";
      fprintf(f09,s500);
      fprintf(f09,"------ ");

      s.flush();
      return;
      }
//=============================
}


//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int  ISUM(int I, int J, int L1, int L2) {

      int ISUM=0;
      for(int L=L1; L<=L2; L++) ISUM+=NR[I][J][L];
      return ISUM;
}


