#include <QDebug>

#include "WinLise_Constant.h"
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define dr13 (1./3.)
#define pow_dr13(A)    pow(A,dr13)
#define pow_mdr13(A)   pow(A,-dr13)
#include "ftype.h"
#include "pace.h"

double RANF();
void PLM();
void FACLOG();
double BASS(int IA1, int IZ1, int IA2, int IZ2, double E, double &Vfusion, double &Rfusion, QTextStream &s);
void TLOM(int IPOT, double AMP, double AZP, double AMT,
          double AZT, double ECM, int &LL, double TL[Max_MOM+1]);
double GetTCN0(double Ecm, double V0, double omega);
double GetEnergyFromL(double Ecm, int L, double reduced_mass, double Rf);
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

extern double MASSES(int IZ, int IN, int &opt);
//extern void COMPND(int IZ, int IN, double ENERGY, int MAXC, QTextStream& s);
extern FILE *f09;
extern const char *s40star;
extern const char *s10star;
extern double MyRandom(void);

extern int _BarrierMode;
extern double _h_omega;

char s35x[10]="         ";
char s100m[77]="----------------------------------------------------------------------------";

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double RANF() {
  // new

  return MyRandom();
  /*
up to 4.18
      const int xran  = 30031961;
      double ranf=random(xran)/double(xran);
//      double ranf = 0.3;

      return ranf;*/
  //  return 0.5;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void PLM() {
  //     WRITTEN BY A. GAVRON
  double PL[25][25];
  //     CALCULATE ASSOCIATED LEGENDRE FUNTIONS AND FORM TABLE OF RUNNING
  //     SUMS FOR LATER SELECTION OF ANGLE.
  //     TO SAVE SPACE THE RUNNING SUM  SUM(N,IANG)  IS STORED FOR EACH
  //     P(L,M) AT  N=L*(L-1)/2+M . THE POLYNOMIAL IS OF ORDER (L-1,M-1).
  //     IANG IS A COS(ANGLE) INDEX FROM .05 TO .95 IN TEN STEPS.
  // ******** WARNING. FOR  M = L  AND L ABOVE 20 RESULTS ARE INACCURATE.
  // ******** REAL*8 SHOULD BE CHANGED TO REAL*16 AND DSQRT BY SQRT IN
  // ******** FORTRAN H EXTENDED, IF ANGULAR DISTRIBUTION AT THESE SPINS
  // ******** IS IMPORTANT.
  extern double _SUM[11][301];//,_CT;      //TRANSPORTED
  // extern int _MDIR;

  double Z,TEMP, PLMM;
  int L,M, LL, MM, N, IANG, IITER;
  double BIN[4]={0,.01667,.05000,.08333};

  FACLOG();
  //     THIS FORMS ALOG(FACTORIALS) FOR LATER  3-J  CALCULATIONS

  for(int IK=1; IK<=300; IK++)
    for(int JK=1; JK<=10; JK++)
      _SUM[JK][IK]=0;

  //==================================================== start cikl 1
  for(IANG=1; IANG<=10; IANG++)
    for(IITER=1; IITER<=3; IITER++) {
        Z=IANG*0.1-BIN[IITER];
        PL[1][1]=1.;
        PL[2][1]=Z;
        TEMP=1.0-pow(Z,2.);
        PL[2][2]=sqrt(TEMP);
        PL[3][1]=(3.*pow(Z,2)-1.0)*0.5;
        PL[3][2]=PL[2][2]*3.*Z;
        PL[3][3]=3.*TEMP;

        for(L=4; L<=24; L++) {
            LL=L-1;
            PL[L][1]=((2*LL-1)*Z*PL[L-1][1]-(LL-1)*PL[L-2][1])/LL;
            PL[L][2]=((2*LL-1)*Z*PL[L-1][2]-LL*PL[L-2][2])/(LL-1);
            for(M=3; M<=L; M++) {
                MM=M-1;
                PL[L][M]=(LL+MM-1)*PL[L-1][M-1]-(LL-MM+1)*Z*PL[L][M-1];
                PL[L][M]=PL[L][M]/PL[2][2];
                //     INDEX OF (L,M) = N
              }
          }

        //          DO 3 L=1,6
        //   3 IF(IITER.EQ.2)WRITE(9, 4)IANG,LL,(PL(L,M),M=1,L)
        //   4 FORMAT(2I10,6F10.4)

        for(L=1; L<=24; L++)
          for(M=1; M<=L; M++) {
              PLMM=PL[L][M]*1.E-16;
              N=L*(L-1)/2+M;
              _SUM[IANG][N]+=pow(PLMM,2);        //TRANSPORTED
            }

      } //continue; L1:
  //==================================================== stop cikl 1

  for(N=1; N<=300; N++) {
      for(IANG=2; IANG<=10; IANG++) _SUM[IANG][N]+=_SUM[IANG-1][N];    //TRANSPORTED
      for(IANG=1; IANG<=10; IANG++) _SUM[IANG][N]/=_SUM[10][N];        //TRANSPORTED
    }

  //          DO 10 L=1,10
  //          DO 10 M=1,L
  //     N=L*(L-1)/2+M
  //  10 WRITE(9, 11)L,M,N,(SUM(N,I),I=1,10)
  //  11 FORMAT(1X,3I4,10F7.4)
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void FACLOG(void) {

  extern double FA[204];

  FA[1]=0.;
  FA[2]=0.;

  for(int I=3; I<=203; I++) FA[I]=FA[I-1]+log(I-1);
}

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
int PACE::COMPOS(int &IZ, int &IN, double &EEXCN, double &VZC, int &MAXC,
                 double &CS_CN, double &ECLOSS, int &INPUT, QTextStream& s)
{
  //
  //     MODIFIED JULIAN SUBROUTINE
  //
  //     ***********
  //=========      common/POT/

  extern double _V0D[5],_R0RD[5],_ARD[5],_R0CD[5]/*,_RCCC[4],_VQS[4],_VQ0[4]*/,_W0D[5],_R0ID[5],_AID[5],
      _RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
  extern int _NPD[5],_IMAG[5],_IRAD[5];

  //=========      common/CP/

  extern double _ALEV[Max_MOM];
  double *SIGJC;
  SIGJC=_ALEV;

  //=========      common/XQSHL/
  //  extern int _NOSHL;
  //=========      common/OLEG/
  // extern int _INPUT, _IDIST;
  extern int _IZP,_IAP,_IZT,_IAT, _IZC, _IAC;
  extern double _SPr, _ST, _QCN;
  extern double _ELAB,_ExpSig,_AGRAZ,_ELOSS ;
  extern int _LMINN, _JCMAX;
  extern double _EEXCN, _EREC, _AJNUC;

  //============

  double TLCN[Max_MOM+1],AJCN[Max_MOM+1];
  double  ELB1, ROT;
  int IA, INP, INT, IDENT, ISPIN, LL, L1=Max_MOM, I, IJUMP,N, N2;
  double AZP, AMP, AZT, AMT, AMC, ECM, EC1,EX, E;
  double SIGOM,PLSQ;
  double FAC,EXPS,TLLS,TLS, Y, LdiffOM=0, LgrazFus=0,SL,SL2;
  int L0=0, L2, ITER,ISC,LMAX,LMIN=0,Li, Lyrast,K,LM1;
  double AL,DEN,LgrazOM=0, RSIG, SF, EXPON, v,SC;
  double SCHX,SCHN,SINB;
  double AC1, AJC, SUMTL;
  double DS;
  double EREC,VZCC,VZB,VZBC;
  double L;
  double Lmin=0;
  double Lmax=Max_MOM;
  double v1,v2;
  double AMAXC,AJNUC;
  double ReducedMass, AKr;
  double Vfusion, Rfusion;

  char format1346[123]="\\par %s Warning! Input channel L above %d"
                       "\\par %s Optical model L-diffuseness ignored."
                       "\\par %s Value set internally to %5.2f";
  char format1346html[132]="<p> %s Warning! Input channel L above %d</p>"
                           "<p>% s Optical model L-diffuseness ignored.</p>"
                           "<p> %s Value set internally to %5.2f</p>";

  if(INPUT>1)goto L1401;
  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

  //   AMP     IS A OF PROJECTILE
  //   AZP     IS Z OF PROJECTILE
  //   SP      IS SPIN OF PROJECTILE
  //   ST      IS SPIN OF TARGET
  //   QCN     IS Q OF COMPOUND NUCLEUS FORMATION
  //   IPOT    IS OPTICAL MODEL POTENTIAL CODE
  //   JCMAX   IS CUT OFF FROM LEVEL DENSITY TABLE
  //
  //     **************************
  AZP=_IZP;        AMP=_IAP;
  AZT=_IZT;        AMT=_IAT;
  INP=_IAP-_IZP;
  INT=_IAT-_IZT;
  IN = INP+ INT;
  _IZC=IZ=_IZP+_IZT;
  _IAC=IA= IZ + IN;
  AMC=IA;

  int opt1, opt2, opt3;
  if(_QCN==0.)  _QCN=MASSES(_IZP,INP,opt1)+MASSES(_IZT,INT,opt2)-MASSES(IZ,IN,opt3);

  IDENT=(_IZP ==_IZT && _IAP == _IAT  ? 1 : 0 );
  ISPIN=(_SPr != 0.  || _ST  !=  0.   ? 1 : 0 );

  if(IDENT==1 && ISPIN==1){
      fprintf(f09,"\\par******* Warning-non zero spins. TLOM needs correction"
                  " for  identical particles here");
      s << "<p>******* Warning-non zero spins. TLOM needs correction"
           " for  identical particles here </p>";
    }
  ELB1=_ELAB-0.5*_ELOSS;

  SC=BASS(_IAP,_IZP,_IAT,_IZT, ELB1, Vfusion, Rfusion, s);
  if(_ExpSig==0.) _ExpSig=max(1e-20,SC);

  //     AGRAZ IS DIFFUSENESS OF TRANSMISSION COEFFICIENTS.
  //     IF AGRAZ.NE.0  TLOM IS BYPASSSED

  //     EXPSIG IS CORRECT EXPERIMENTAL XSECTION IF EXISTS. IT IS USED
  //     TO CUT THE TL DISTRIBUTION TO PROVIDE THE CORRECT XSECTION.
  //     JCMAX IS MAXIMUM J OF C.N. (FOR CUTOFF)

  ECM=_ELAB*AMT/(AMP+AMT);
  ECLOSS*=AMT/(AMP+AMT);
  EC1=ECM-ECLOSS*0.5;
  EEXCN=_QCN+ECM;
  ReducedMass = AMP*AMT/(AMP+AMT);
  AKr=2.*aem_MeV * ReducedMass*EC1 / (hc_MeV_fermi*hc_MeV_fermi);
  PLSQ=(AKr==0 ? 0 : 10.0*PI/AKr );

  if(EEXCN<=0) {
      fprintf(f09,"\\par\\par\\b  Excitation energy (%.2f) less than 0!\\b0 ",EEXCN);
      s << "<p>&nbsp;</p><p><b> Excitation energy (" << QString::number(EEXCN) << ") less than 0!</p>";
      return -4;
    }
  ROT=34.548*pow(AMC,(-1.6667));
  Lyrast=sqrt(EEXCN/ROT)+.5;  // Lyrast IS SPIN OF YRAST (RIGID BODY) LEVEL AT C.N. EXCITATION  ENERGY.

  // ==================================== Quantum-Mech   START
  if(_BarrierMode==777) { //==============================           was to 1    -- modified 5/22/2013
      //                                          in order to pass Classic "Momentum distribution"
SpecialModeForQM:

      LL=Max_MOM;
      L1=Max_MOM;

      SIGOM=0.;


      for(I=1; I<=Max_MOM; I++) {

          if(IDENT!=0) FAC=(I%2 == 1 ? 2 : 0);
          else         FAC=1;

          if(FAC==0) { TLCN[I]=0; continue;}

          E = GetEnergyFromL(EC1, I-1, ReducedMass, Rfusion);
          TLCN[I]= GetTCN0(E, Vfusion, _h_omega);

          SIGOM+=(2.*I - 1.)*TLCN[I]*FAC;
        }

      SIGOM*=PLSQ;

      goto L7;
    }
  // ==================================== Quantum-Mech  END




L_REPEAT:
  //================================================================ _AGRAZ  start

  if(_AGRAZ!=0.) {
      LL=Max_MOM;
      L1=Max_MOM;

      for(I=1; I<=Max_MOM; I++) {

          EX = (I-40.)/_AGRAZ;         // why 40 is taken?   OT

          if(EX >  20.) TLCN[I]= 0.;
          else if(EX < -20.) TLCN[I]= 1.;
          else               TLCN[I]= 1./(1.+exp(EX));

          //     SADDLE POINT LEVEL DENSITY NOW RIGOROUSLY CALCULATED.
          //     PREVIOUSLY, G.S. LEVEL DENSITY WAS EXPONENTIATED BY SQRT(AF/AN)
          //     FISSION PROBABILITIES MAY DIFFER BY 20 PERCENT.
        }
    }
  else {
      TLOM(4,AMP,AZP,AMT,AZT,EC1,LL,TLCN);
      if(_JCMAX==0)_JCMAX=LL;
      L1=min(LL,_JCMAX)+1;
      if(L1>Max_MOM) {
          fprintf(f09,"\\par  Dimension AF AJCN too small - %i required",L1);
          s << "<p>   Dimension AF AJCN too small - " << QString::number(L1) << " required</p>";
          s.flush();
          return -3;
        }
      for( I=1; I<=L1; I++)
        if(TLCN[I]<0) {
            _AGRAZ=4;
            fprintf(f09,"\\par\\b       Error in Optical model! Probably excitation energy is high ");
            fprintf(f09,"\\par       AGRAZ(diffusness) is used equal to 4 \\b0 ");
            s << "<b><p>      Error in Optical model! Probably excitation energy is high </p><p>     AGRAZ(diffusness) is used equal to 4 </p></b>";
            s.flush();
            goto L_REPEAT;
          }
    }

  //================================================================ _AGRAZ  stop

  //================================================================= SIGOM start
  SIGOM=0.;

  for( I=1; I<=L1; I++) {
      if(IDENT!=0) FAC=(I%2 == 1 ? 2 : 0);
      else         FAC=1;

      SIGOM+=(2.*I - 1.)*TLCN[I]*FAC;
    }

  SIGOM*=PLSQ;
  //================================================================= SIGOM stop

  //================================================================= EXPSIG!=0  start
  if(_ExpSig==0.)goto L7;
  //
  //   NORMALIZATION OF OPTICAL MODEL AND EXPERIMENTAL CROSS SECTIONS
  //
  //   T = 1/(1+EXP((L-L(G))/V)
  //   M*L+B = Y WHERE Y = LN(1/T-1), M = 1/V, AND B = -M * L(G)
  //         SOLVE FOR L(G) AND V FROM CALCULATED TRANSMISSION COEFFICIENTS
  //   S = C * SUM((2 * L + 1) * T  WHERE L(G) IN T IS NOW L(C)
  //   (L(G)+1)**2/(L(C)+1)**2 = SIGMA(CALCULATED)/SIGMA(EXPERIMENTAL)
  //         IS USED TO SOLVE FOR L(C)
  //

  EXPS=_ExpSig/PLSQ;
  TLLS=0.;
  TLS=0.;
  IJUMP=0;

  L2=L1;

  for(I=2; I<=L1; I++) {      //-----------------------------------------

      L=I-1;
      if(TLCN[I]> 0.9*TLCN[1]) continue;
      if(TLCN[I]< 0.004) {
          L2=L;
          break;
        }

      if(IJUMP!=1) {
          IJUMP=1;
          L0=L;
        }

      Y=log(1./TLCN[I]-1.);
      TLS  +=Y;
      TLLS +=Y*L;
    }                                      //  L0 and L2 boundaries
  //-----------------------------------------


  //-----------------------------------------
  if(TLLS<=0.) {
      fprintf(f09,"\\par %s Max. spin probably too small (MAXC,AJNUC) "
                  "\\par %s  or reaction cross-section too big . Change AGRAZ",s10star,s10star);
      s << "<p> " << QString::fromUtf8(s10star) << "Max. spin probably too small (MAXC,AJNUC)  </p><p> or reaction cross-section too big . Change AGRAZ" << QString::fromUtf8(s10star) << "</p>";
      _AGRAZ=0.5;
      LdiffOM=75.;
      LgrazFus=_AGRAZ;
      fprintf(f09,format1346,s10star,Max_MOM,s10star,s10star,_AGRAZ);
    }
  else {
      double L2_Sq = (L2-1)*L2/2;       double   L2S =  L2_Sq*(2*L2-1)/3;
      double L0_Sq = (L0-1)*L0/2;       double   L0S =  L0_Sq*(2*L0-1)/3;
      L2_Sq -= L0_Sq;                     L2S -= L0S;

      SL   = L2_Sq;
      SL2 = L2S;
      AL  = L2-L0;

      DEN = SL*SL-AL*SL2;
      v1  = TLS*SL  - AL*TLLS;       if(fabs(v1)<1e-10)v1=1;
      v2  = TLS*SL2 - SL*TLLS;       if(fabs(v2)<1e-10)v2=1;

      LdiffOM  = DEN / v1;       ///  LdiffOM == ELAM

      LgrazFus = v2/v1;      // LgrazFus == ELG

    }
  //-----------------------------------------



  LgrazOM=LgrazFus;
  ITER=0;

  //=====================================

  for(ITER=1; ITER<=100; ITER++)
    {
      SF=0;

      for(I=1; I<=Max_MOM; I++)
        {

          L=I-1;
          EXPON=(L-LgrazFus )/LdiffOM;
          if(EXPON>40.)   break;   // before it was break

          DEN=1.+ (EXPON <= -40. ? 0 : exp(EXPON));
          v=(2*L+1)/DEN;
          if(  v < 0.001)  break;
          SF+=v;
        }

      RSIG=EXPS/SF;
      if(fabs(RSIG-1.)<0.0001) break;

      if(RSIG>1)  Lmin=LgrazFus;
      else        Lmax=LgrazFus;

      LgrazFus=(Lmax+Lmin)/2.;
    }

  //=============================================================
  if(ITER==101 && RSIG < 1) {

      // NEW   05/22/2013  ------------------------------------------------
      if(_BarrierMode==1) goto SpecialModeForQM;

      fprintf(f09,"\\par******* \\b Diffusness has been changed. Cross-section value is vey low! (OT)\\b0 ");
      fprintf(f09,"\\par******* \\b Try to use quantum-mechanical model! (OT)\\b0 ");
      s << "<p><b> *******  Diffusness has been changed. Cross-section value is vey low! (OT)</b></p><p><b>*******  Try to use quantum-mechanical model! (OT)</b></p>";
      Lmin=0;
      LdiffOM*=2;
      Lmax=LdiffOM;

      for(ITER=1; ITER<=100; ITER++)  {
          SF=0;
          for(I=1; I<=Max_MOM; I++) {

              L=I-1;
              EXPON=(L-LgrazFus )/LdiffOM;
              if(EXPON>40.)   break;   // before it was break

              DEN=1.+ (EXPON <= -40. ? 0 : exp(EXPON));
              if( (v=(2*L+1)/DEN) < 0.001)  break;
              SF+=v;
            }

          RSIG=EXPS/SF;
          if(fabs(RSIG-1.)<0.0001) break;

          if(RSIG>1)  Lmin=LdiffOM;
          else        Lmax=LdiffOM;

          LdiffOM=(Lmax+Lmin)/2.;
        }
    }
  //=============================================================
  if(ITER==101 && RSIG < 1){
      fprintf(f09,"\\par******* \\b Transmission coeficients were normalized. Cross-section value is vey low! (OT)\\b0 ");
      s << "<p><b>*******  Transmission coeficients were normalized. Cross-section value is vey low! (OT)</b></p>";
    }else RSIG=1;


  //=============================================================
  for( I=1; I<=Max_MOM; I++) {
      L=I-1;
      EXPON = min(60,max(-60,(L-LgrazFus)/LdiffOM)); //  60 DEPENDS ON COMPUTER IN USE.
      if(  (TLCN[I]=RSIG/(1.+exp(EXPON))) < 0.001) { L1=L;  break; }
    }
  //================================================================= EXPSIG!=0  stop


  //=============================================================

L7:   CS_CN=0;


  if( !(L1<Max_MOM || _AGRAZ!=0.)) {
      _AGRAZ=0.5;
      fprintf(f09,format1346,s10star,s10star,s10star,_AGRAZ);
      s << format1346html;
      L1=Max_MOM;
    }
  SCHX=_SPr+_ST;
  SCHN=fabs(_SPr-_ST);
  SINB=(2.0*_SPr+1.0)*(2.0*_ST+1.0);
  ISC=_SPr+_ST+0.01;
  SC=_SPr+_ST - ISC;

  //---------------------------------------

  for(I=1; I<=L1; I++) {//41
      AC1=I-1;
      AJC=AC1+SC;
      SUMTL=0.;
      LMAX=(AJC+SCHX)+0.01;

      if(SCHX < AJC)    LMIN = AJC  - SCHX + 0.01;
      else if(SCHN > AJC)    LMIN = SCHN - AJC  + 0.01;
      else if(SCHN ==AJC)    LMIN = 0;

      Li=LMIN;

      while( !(Li>LL||Li>L1)  ) {

          DS= min(SCHX,     AJC+Li )
              -max(SCHN,fabs(AJC-Li))
              +1;

          //        if(DS >= 1.)  SUMTL += DS*TLCN[Li+1];
          if(DS >= 1. && Li<L1)  SUMTL += DS*TLCN[Li+1];    //  Oleg 19 may 2003
          if(++Li>LMAX) break;
        }

      AJCN[I]=AJC;
      SIGJC[I]=(2.0*AJC+1.0)*PLSQ*SUMTL/SINB;
      if(IDENT!=0){
          if(I%2==1)SIGJC[I]*=2.;
          else      SIGJC[I] =0.;
        }
      //      PARTIAL XSECTION (IS  ALEV  IN MAIN PROGRAM)
      CS_CN+=SIGJC[I];
    }
  //---------------------------------------

  EREC=_ELAB*AMP/AMC;
  VZC=sqrt(2.0*EREC/(aem_MeV*AMC));
  VZCC=VZC*30.;
  VZB=sqrt(2.0*_ELAB/(aem_MeV*AMP));
  VZBC=VZB*30.;
  MAXC=L1;

  fprintf(f09,"\\f0\\fs20\\par\\par%s",s100m);
  fprintf(f09,"\\par\\cf3\\b      Starting conditions\\cf0\\b0");
  s << "<h3 style=\"color:blue\">    Starting conditions </h3><table cellpadding=15><tr>";
  fprintf(f09,"\\par%s                     Z     N     A    Spin", s35x);
  s << "<th></th><th>Z</th><th> N</th><th>A</th> <th>Spin</th></tr>";
  fprintf(f09,"\\par%s\\i Projectile\\i0          %3d   %3d   %3d   %5.1f",s35x,_IZP, INP, _IAP,_SPr);
  s << "<tr><td><em>Projectile</em></td><td align=center>" << QString::number(_IZP) << " </td><td align=center>" << QString::number(INP) << " </td><td align=center>" << QString::number(_IAP) << "</td><td align=center>" << QString::number(_SPr) << "</td></tr>";
  fprintf(f09,"\\par%s\\i Target    \\i0          %3d   %3d   %3d   %5.1f",s35x,_IZT, INT, _IAT,_ST);
  s << "<tr><td><em>Target</em></td><td align=center>" << QString::number(_IZT) << "</td><td align=center>" << QString::number(INT) << "</td><td align=center>" << QString::number(_IAT) << "</td><td align=center>" << QString::number(_ST) << "</td></tr>";
  fprintf(f09,"\\par%s\\i Compound nucleus\\i0    %3d   %3d   %3d \\par\\fs16",        s35x, IZ , IN , IA);
  s << "<tr><td><em>Compound nucleus</em></td><td align=center>" << QString::number(IZ) << " </td><td align=center>" << QString::number(IN) << " </td><td align=center> " << QString::number(IA) <<  "</td></tr></table>";
  fprintf(f09,"\\par \\i  Bombarding energy (MeV)                    \\i0  %10.2f",_ELAB);
  s << "<table><tr><td><em> Bombarding energy (MeV)</em></td><td>  " << QString::number(_ELAB,'f',2) << " </td></tr>";
  fprintf(f09,"\\par \\i  Center of mass energy (MeV)                \\i0  %10.2f",ECM);
  s << "<tr><td><em> Center of mass energy (MeV)</em> </td><td>" << QString::number(ECM,'f',2) << "</td></tr>";
  fprintf(f09,"\\par \\i  Compound nucleus excitation energy (MeV)   \\i0  %10.3f",EEXCN);
  s << "<tr><td><em>   Compound nucleus excitation energy (MeV)  </em></td><td>  " << QString::number(EEXCN,'f',3) << "</td></tr>";
  fprintf(f09,"\\par \\i  Q-value of reaction (MeV)                  \\i0  %10.3f",_QCN);
  s << "<tr><td><em>   Q-value of reaction (MeV)   </em></td><td>  " << QString::number(_QCN,'f',3) << "</td></tr>";
  if(ECLOSS>0.){
      fprintf(f09,"\\par \\i  Excitation energy loss thru target (MeV)\\i0     %10.2f",ECLOSS);
      s << "<tr><td><em> Excitation energy loss thru target (MeV) </em></td><td> " << QString::number(ECLOSS,'f',3) << "</td></tr>";
    }
  fprintf(f09,"\\par \\i  Compound nucleus recoil energy (MeV)     \\i0    %10.3f",EREC);
  s << "<tr><td><em>   Compound nucleus recoil energy (MeV)    </em></td><td> " << QString::number(EREC,'f',3) << "</td></tr>";
  fprintf(f09,"\\par \\i  Compound nucleus recoil velocity (cm/ns) \\i0    %10.3e", VZCC);
  s << "<tr><td><em>   Compound nucleus recoil velocity (cm/ns) </em></td><td>  " << QString::number(VZCC,'e',3) << "</td></tr>";
  fprintf(f09,"\\par \\i  Compound nucleus velocity/c \\i0                 %10.3e", VZC);
  s << "<tr><td><em> Compound nucleus velocity/c </em> </td><td> " << QString::number(VZC,'e',3) << "</td></tr>";
  fprintf(f09,"\\par \\i  Beam velocity (cm/ns) \\i0                       %10.3e",VZBC);
  s << "<tr><td><em>  Beam velocity (cm/ns)  </em> </td><td> " << QString::number(VZBC,'e',3) << "</td></tr>";
  fprintf(f09,"\\par \\i  Beam velocity/c    \\i0                          %10.3e",VZB);
  s << "<tr><td><em>Beam velocity/c  </em></td><td>  " << QString::number(VZB,'e',3) << "</td></tr></table>";
  if(_BarrierMode==0)
    {
      if(_AGRAZ==0.) {
          fprintf(f09,"\\fs14\\par\\par\\b  Optical model parameters for entrance channel\\b0"
                      "\\par   V       *E    *E**2   R0R   ARD    R0C     W0    *E   *E**2  R01    AID    RMCHD    NPD    IRAD"
                      "\\par %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7d %4d %4d",
                  _V0D[4],_V01D[4],_V02D[4],_R0RD[4],_ARD[4],_R0CD[4],_W0D[4],_W01D[4],_W02D[4],_R0ID[4],_AID[4],
              _RMCHD[4],_NPD[4], _IMAG[4],_IRAD[4]);
          s << "<p>&nbsp;</p><table cellpadding='5'><tr><th colspan=15 align=center> Optical model parameters for entrance channel </th></tr>"
               "<tr><th>V</th><th>*E</th><th>*E**2</th><th>R0R</th><th>ARD</th><th>R0C</th><th>W0</th><th>*E</th><th>*E**2</th><th>R01</th><th>AID</th><th>RMCHD</th><th>NPD</th><th>IRAD</th></tr>";
          s << "<tr><td>" << QString::number(_V0D[4]) << "</td><td>" << QString::number(_V01D[4]) << "</td><td>" << QString::number(_V02D[4]) << "</td><td>"  <<
                                                                                                                                                 QString::number(_R0RD[4]) << "</td><td>" << QString::number(_ARD[4]) << "</td><td>" << QString::number(_R0CD[4]) << "</td><td>" << QString::number(_W0D[4]) << "</td><td>"
                                                                                                                                                                                                                                                                                                             << QString::number(_W01D[4]) << "</td><td>" << QString::number(_W02D[4]) << "</td><td>" << QString::number(_R0ID[4]) << "</td><td>" << QString::number(_AID[4])
              << "</td><td>" << QString::number(_RMCHD[4]) << "</td><td>" << QString::number(_NPD[4]) << "</td><td>" << QString::number(_IMAG[4]) << "</td><td>" << QString::number(_IRAD[4])
              << "</td></tr></table>";

          fprintf(f09,"\\par"
                      "\\par\\fs16\\i  Optical model reaction cross section (mB)\\i0  %10.3f",SIGOM);
          s << "<p>&nbsp;</p><p><em> Optical model reaction cross section (mB)</em> " << QString::number(SIGOM) << "</p>";
        }
      else   {
          fprintf(f09,"\\par\\cf2\\fs20\\f2"
                      "\\par *** Input transmission coefficients determined by input value of TL diffuseness."
                      "\\par *** diffuseness = %8.2f"
                      "\\par *** Optical model input calculation bypasses. *********\\fs16\\cf0\\f0\\par",_AGRAZ);
          s << "<p>&nbsp;</p><p> *** Input transmission coefficients determined by input value of TL diffuseness. </p>"
               "<p> *** diffuseness = " << QString::number(_AGRAZ,'f',2) << "<p>*** Optical model input calculation bypasses. *********</p> ";
        }
    }
  s << "<p>&nbsp;</p>";
  if(_ExpSig!=0.) {
      s << "<table>";
      if(_AGRAZ==0.) {
          fprintf(f09,"\\par%s\\i OM L-grazing  \\i0     %10.3f  "   ,s35x,LgrazOM);
          s << "<tr><td><em> OM L-grazing </em></td><td> " << QString::number(LgrazOM) << " </td></tr>";
          fprintf(f09,"\\par%s\\i OM L-diffuseness \\i0  %10.3f\\par",s35x,LdiffOM);
          s << "<tr><td><em> OM L-diffuseness </em></td><td> " << QString::number(LdiffOM) << "</td> </tr>";
        }
      if(_LMINN>0){
          fprintf(f09,"\\par \\i Partial waves below L =\\i0 %-3d\\i exluded\\i0",_LMINN);
          s << "<tr><td><em>  Partial waves below L = </em></td><td>" << QString::number(_LMINN) << " excluded</td></tr>";
        }
      fprintf(f09,"\\par\\par \\i  Experimental fusion cross section (mb) \\i0      %10.3g",_ExpSig);
      s << "<tr><td><em> Experimental fusion cross section (mb)</em> </td><td>  " << QString::number(_ExpSig,'e',2) << "</td></tr>";
      if(_BarrierMode==0) {
          fprintf(f09,"\\par \\i  Fusion L-grazing                       \\i0      %10.3f",LgrazFus);
          s << "<tr><td><em>  Fusion L-grazing </em></td><td>" << QString::number(LgrazFus,'f',2) << "</td></tr>";
          fprintf(f09,"\\par \\i  Fusion L-difusseness                   \\i0      %10.3f",LdiffOM);
          s << "<tr><td><em>  Fusion L-diffuseness </em></td><td>" << QString::number(LdiffOM,'f',2) << "</td></tr>";
        }
      fprintf(f09,"\\par \\i  Yrast spin at maximum excitation energy \\i0     %10d"  ,Lyrast);
      s << "<tr><td><em>  Yrast spin at maximum excitation energy   </em></td><td>" << QString::number(Lyrast) << "</td></tr>";
      fprintf(f09,"\\par \\i  Compound nucleus formation cross section (mb)\\i0 %10.3g",CS_CN);
      s << "<tr><td><em>   Compound nucleus formation cross section (mb) </em></td><td> " << QString::number(CS_CN,'e',2) << "</td></tr>";
    }

  s.flush();
  SC=((IA)%(2))*0.5;

  if(_LMINN==0){s << "</table>"; goto L45;}



  for(LL=1; LL<=_LMINN; LL++) {  //    ANG. MOM. LMINN LESS 1 IS INDEX LMINN
      CS_CN-=SIGJC[LL];
      SIGJC[LL]=0.;
    }
  if (_ExpSig==0) {s << "<table>";}
  fprintf(f09,"\\par  Compound nucleus cross section in window (mb) %10.3f",CS_CN);
  s << "<tr><td><em>  Compound nucleus cross section in window (mb) </em></td><td>" << QString::number(CS_CN) << "</td></tr></table>";
L45:  fprintf(f09,"\\par\\par%s",s100m);
  //s << "<p>&nbsp;</p><hr />";// << QString::fromUtf8(s100m) << "</p>";
  fprintf(f09,"\\par%s\\b           Partial cross sections (mb)\\b0", s35x);
  //  s << "<h3>  Partial cross sections (mb) </h3>";
  fprintf(f09,"\\par%s",s100m);
  //s << "<p>" << QString::fromUtf8(s100m) << "</p>";
#define COL  5

  fprintf(f09,"\\par ");
  s << "<table border=1 cellspacing='0' cellpadding='5'><tr><th colspan=10>Partial cross sections (mb)</th></tr><tr>";
  for(K=1; K<=COL; K++) {
      fprintf(f09,"   J   SIG(J) |");
      s << "<th align=center>   J </th><th align=center>  SIG(J)  </th>";
    }
  fprintf(f09,"\\par%s",s100m);
  s <<"</tr>";
  N2=L1/COL+ (L1%COL > 0 ? 1 : 0 );

  for(I=1; I<=N2; I++) {
      fprintf(f09,"\\par ");
      s << "<tr>";
      for(K=1; K<=COL; K++) {
          N=I+N2*(K-1);
          if(N<=L1) {
              fprintf(f09,"\\cf3 %5.1f\\cf0  %7.2g |",AJCN[N],SIGJC[N]);
              s << "<td align=center style=\"color:blue\">" << QString::number(AJCN[N]) << "</td><td align=center> " << QString::number(SIGJC[N],'g',2) << "</td>";
            }else {
              fprintf(f09,"              |");
              s <<  "<td></td><td></td>";
            }
        }
      s << "</tr>";
    }
  fprintf(f09,"\\par%s\\par",s100m);
  s << "</table><p>&nbsp;</p>";// << QString::fromUtf8(s100m) << "</p>";
  return 1;


  //WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  //------------------------------   INPUT > 1

L1401:

  //      AJNUC IS REAL TRUE SPIN. (NOT THE INDEX)
  // 1401  READ(10,1402)IZ,IA,EEXCN,EREC,AJNUC ,LMINN
  //      type *,IZ,IA,EEXCN,EREC,AJNUC ,LMINN

  IZ=_IZC; IA=_IAC; // Oleg
  EEXCN = _EEXCN;   // Oleg
  EREC = _EREC;     // Oleg
  AJNUC = _AJNUC;   // Oleg

  IN= IA-IZ;
  AMC=IA;
  VZC=sqrt(2.0*EREC/(aem_MeV*AMC));
  VZCC=VZC*30.;
  SC=((IA)%(2))*0.5;
  AMAXC=AJNUC -SC+1.1;
  MAXC=AMAXC;

  //O//LEG      if(INPUT==3)fscanf (10,1403)MAXC,(SIGJC[I],I=1,MAXC);
  if(INPUT==4)COMPND(IZ,IN,EEXCN,MAXC, s);

  //L1403: format(I5,5X,7F10.5/(8F10.5));
  fprintf(f09,"\\fs20\\par\\par%s",s100m);
  s << "<p>&nbsp;</p><p>";
  fprintf(f09,"\\par\\b%s  Starting compound nucleus\\b0", s35x);
  s << "<hr /><h3 style=\"color:blue\">  Starting compound nucleus </h3> ";// << QString::fromUtf8(s35x) << "</p>";
  //      fprintf(f09,"\\par%s",s100m);


  fprintf(f09,"\\par%s%s \\b                     Z     N     A\\b0", s35x, s35x);
  s << "<table><tr><th></th><th>Z</th><th>N</th><th>A</th></tr>";
  fprintf(f09,"\\cf3\\par%s%s\\i\\b Compound nucleus\\i0     %3d   %3d   %3d \\cf0\\b0\\fs16\\par", s35x, s35x, IZ , IN , IA);
  s << "<tr><td><em> Compound nucleus </em></td><td> " <<  QString::number(IZ)
    << "</td><td> " << QString::number(IN) << "</td><td> " << QString::number(IA) << "</td></tr></table><p>&nbsp;</p>";
  fprintf(f09,"\\par\\i  Compound nucleus excitation energy (MeV)\\i0   %10.2f",EEXCN);
  s << "<table><tr><td><em> Compound nucleus excitation energy (MeV) </em></td><td>" << QString::number(EEXCN) << "</td></tr>";
  fprintf(f09,"\\par\\i  Compound nucleus recoil energy (MeV)\\i0       %10.2f",EREC);
  s << "<tr><td><em> Compound nucleus recoil energy (MeV) </em></td><td>" << QString::number(EREC) << "</td></tr>";
  fprintf(f09,"\\par\\i  Compound nucleus recoil velocity (cm/ns)\\i0   %10.2e", VZCC);
  s << "<tr><td><em> Compound nucleus recoil velocity (cm/ns) </em></td><td>" << QString::number(VZCC) << "</td></tr>";
  fprintf(f09,"\\par\\i  Compound nucleus velocity/c\\i0                %10.2e", VZC);
  s << "<tr><td><em>  Compound nucleus velocity </em></td><td>" << QString::number(VZC) << "</td></tr>";

  if(INPUT<3) {
      fprintf(f09,"\\par\\i  Compound nucleus spin\\i0                      %10.1f",AJNUC);
      s << "<tr><td><em>  Compound nucleus spin </em></td><td>" << QString::number(AJNUC) << "</td></tr></table><hr />";
      CS_CN=100.;
      fprintf(f09,"\\fs20\\par\\fc0%s\\fs16",s100m);
      //s << "<p>" << QString::fromUtf8(s100m) << "</p>";
      s.flush();
      return 1;
    }
  s << "</table><hr />";
  CS_CN=0.;
  LM1=_LMINN+1;

  for( I=1; I<=MAXC; I++) {
      AC1=I-1;
      AJCN[I]=AC1+SC;
      if(INPUT==5) {
          if(I>=LM1)SIGJC[I]=2.*AJCN[I]+1.;
          else SIGJC[I]=0.;
        }
      CS_CN+=SIGJC[I];
    }


  fprintf(f09,"\\par\\par%s",s100m);
  // s << "<p>" << QString::fromUtf8(s100m) << "</p>";
  fprintf(f09,"\\par%s\\b           Partial cross sections (mb)\\b0", s35x);
  s << "<h3 align=center> Partial cross sections (mb)</h3><p>&nbsp;</p>";
  fprintf(f09,"\\par%s",s100m);
  // s << "<p>" << QString::fromUtf8(s100m) << "</p><p>";
  s << "<table align=center cellpadding='5'><tr>";
  fprintf(f09,"\\par");
  for(K=1; K<=COL; K++) {
      s << "<td align=center>J</td><td align=center>SIG(J)</td>";
      fprintf(f09,"   J   SIG(J) |");
    }
  s << "</tr>";
  N2=MAXC/COL+ (MAXC%COL > 0 ? 1 : 0 );

  L1=min(L1,MAXC);  // Oleg

  for(I=1; I<=N2; I++) {
      fprintf(f09,"\\par");
      s << "<tr>";
      for(K=1; K<=COL; K++) {
          N=I+N2*(K-1);
          if(N<=L1) {
              s << "<td align=center style=\"color:blue\">" << QString::number(AJCN[N],'f',1) << "</td><td align=center>" << QString::number(SIGJC[N],'g',3) << "</td>";
              fprintf(f09,"\\cf3 %5.1f\\cf0  %7.2g |",AJCN[N],SIGJC[N]);
            }else {
              s << "<td></td>";
              fprintf(f09,"              |");
            }
        }
      s << "</tr>";
    }
  s << "</table>";
  fprintf(f09,"\\par%s\\par",s100m);
  s << "<p>  Compound nucleus formation cross section (mb) " << QString::number(CS_CN) << "</p><p>&nbsp;</p>";
  fprintf(f09,"\\par  Compound nucleus formation cross section (mb) %10.3f",CS_CN);
  s.flush();
  return 1;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
double BASS(int IA1, int IZ1, int IA2, int IZ2, double E, double &Vfusion, double &Rfusion, QTextStream& s)

{
  //     WRITTEN BY A. GAVRON
  //  changed by Oleg

  fprintf(f09,"\\par\\fs18\\i  ************ Fusion cross-section taken from Bass model\\fs16\\i0");
  s << "<p><em> ************ Fusion xsection taken from Bass model</em></p>";
  double Z1=IZ1;
  double Z2=IZ2;
  double A1=IA1;
  double A2=IA2;
  double R10=pow(A1,.3333);
  double R20=pow(A2,.3333);

  double R1=R10*1.16 - 1.39/R10;
  double R2=R20*1.16 - 1.39/R20;

  const double OlegCorrection[24] = {0.00,0.76,0.43,0.29,0.22,0.17,0.13,0.11,0.09,0.08,
                                     0.06,0.05,0.04,0.04,0.03,0.03,0.02,0.02,0.02,0.01,
                                     0.01,0.01,0.01,0.01};                      // difference between Cp and Bass

  if(IA1<=23 && IA1 > 0) R1 += OlegCorrection[IA1];
  if(IA2<=23 && IA2 > 0) R2 += OlegCorrection[IA2];


  double ECM=E*A2/(A1+A2);
  double SIGmin=1e20;
  double S, G, VNUC, VCOUL, V, SIGF, R;
  //      double SIG_Bmax=0;
  Vfusion=0;

  int I0=20.*(R1+R2);

  // looking for Barrier max

  /********** MICHELLE DEBUG for I0 divide by zero error **************/
  if(I0 < 0){

      I0 = 0;
    }
  for(int I=I0; I<=599; I++) {  //L2:
      R=0.05*I;
      S=R-R1-R2;
      if(S<=0.) continue;

      G=1./(0.03*exp(S/3.30)+.0061*exp(S/0.65));
      VNUC=-R1*R2*G/(R1+R2);
      VCOUL=1.44*Z1*Z2/R;
      V=VNUC+VCOUL;
      double CScoef =  R*R*10.;

      if(_BarrierMode==0) SIGF=PI* CScoef* (1.-V/ECM);
      else {
          double t1 = 2.* PI / _h_omega * (ECM-V);
          double t3 = ( t1 > 100 ? t1 :  log(1.+ exp(t1)));
          double t4 = CScoef* _h_omega / 2./ ECM;
          SIGF = t4*t3;
        }

      if( SIGmin > SIGF) {
          Rfusion=R;
          SIGmin=SIGF;
        }

      Vfusion=max(V,Vfusion);
    }
  //============================================ end cikl L2

  if(SIGmin<0.)SIGmin=0.;
  if(_BarrierMode==0 && Vfusion > ECM) SIGF=0;    // Oleg

  fprintf(f09,"\\par\\fs18\\i  Bass fusion xsection for E =\\i0 %6.1f\\i  MeV  is\\i0  %.4g\\i  mb"
              "\\par  Fusion radius =\\i0 %6.2f\\i  fm. Barrier height is\\i0 %7.2f\\i  MeV\\i0"
              "\\par  Transmission probability for a one-dimens.barrier: \\b %s \\b0 \\fs16",
          E,SIGmin,Rfusion,Vfusion,(_BarrierMode ==0 ? "Classical" : "Quantum-Mechanical"));
  s << "<p><em> Bass fusion xsection for E = </em>" << QString::number(E) << "<em> MeV is </em>" << QString::number(SIGmin)
    << "<em> mb </em> </p><p><em> Fusion radius = </em>" << QString::number(Rfusion) << "<em> fm. Barrier height is </em>"
    << QString::number(Vfusion) << "<em> MeV </em></p><p> Transmission probability for a one-dimens. barrier: <b>" <<
       (_BarrierMode ==0 ? "Classical" : "Quantum-Mechanical") << "</b></p>";
  if(_BarrierMode==1){
      fprintf(f09,"\\par\\fs18  h_omega (curvature parameter) = %.1f\\i  MeV\\i0 \\fs16", _h_omega);
      s << "<p>   h_omega (curvature parameter) = " << QString::number(_h_omega) << " <em>MeV</em></p>";
    }
  s << "<hr />";
  s.flush();
  return SIGmin;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//
//     ORIGINAL JULIAN SUBROUTINE
//
//     ***********

void TLOM(int IPOT, double AMP, double AZP, double AMT,
          double AZT, double ECM, int &LL, double TL[Max_MOM+1])
{
  //=========      common/POT/
  extern double _V0D[5],_R0RD[5],_ARD[5],_R0CD[5],_RCCC[4],/*_VQS[4],_VQ0[4],*/_W0D[5],_R0ID[5],_AID[5],
      _RMCHD[5],_V01D[5],_V02D[5],_W01D[5],_W02D[5];
  extern int _NPD[5],_IMAG[5],_IRAD[5];

  //====================================       common/SRCH/ - OLEG
  double VR[501],VI[501],URE[501],UIM[501],FC[Max_MOM+1],GC[Max_MOM+1],
      COE[7],COAA[7],COAC[7],COAY[7];
  //====================================

  double SIG[Max_MOM+1],SIGL[Max_MOM+1];
  double V0, R0R, AR, R0C,W0, R0I, AI, AP3, AT3, RR, RC, RI, RMCH=0;
  double ANP=0, DR=0,ZZO=0,ZZI=0,R=0, EXR=0, EXI=0, VV0=0, VC=0, AM=0, VW0=0; //MPK Initialized these
  int I, IRADUS, IMG,NPINT=0,NP=0,LLX,LX,L2, L1, LL1, L3, L4,N, M, I1, I2;
  double AK,PLSQ,EP,ETAC,Q,DX,ALX,FC0, AL2, FL1, ZZ, AL1;
  double ETA,ETASQ,ETA2,ETA6, RX,A1, A2,B1,B2,FP=0,GP=0,D=0, STEP=0,FL=0;
  double SMAX,COH,COA,COC,COAS,COAB, COAE, COAT, COAA1,COAD,STP,FQ;
  // int    IFN,MM,J, L5,Li, II, MI, K1; // MPK checking memory leaks from Dr. Memory
  int   J, L5,Li, II, MI, K1;
  int IFN = 0, MM = 0;
  double RENORM,FSUBL,GSUBL;
  double GLP,FLP,TEMPA,TEMPB, SSUBL,DELSL,RXX,XII, XX,DX2, DRE, DI;
  double DEN,FLR,FLI, DDD,ETARE, ETAIM, SIGR;
  //   I=IPOT;

  for(int i=0; i<=Max_MOM; i++) GC[i]=0; // Oleg 09/07/22

  if(IPOT==0)goto L10001;

  I=IPOT;
  IRADUS=_IRAD[I];
  V0=_V0D[I]+_V01D[I]*ECM+_V02D[I]*ECM*ECM;
  R0R=_R0RD[I];
  AR =_ARD[I];
  R0C=_R0CD[I];
  if(R0C==0.)R0C=R0R;
  IMG=_IMAG[I];
  W0=_W0D[I]+_W01D[I]*ECM+_W02D[I]*ECM*ECM;
  R0I=_R0ID[I];
  AI =_AID[I];
  AP3=pow(AMP,0.3333);
  AT3=pow(AMT,0.3333);

  switch(IRADUS) {
    case 1 :       RR=R0R*AT3;
      RC=R0C*AT3;
      RI=R0I*AT3;
      break;
    case 2 :        RR=R0R*(AP3+AT3);
      RC=R0C*(AP3+AT3);
      RI=R0I*(AP3+AT3);
      break;
    default :       RR=R0R;
      RC=R0C;
      RI=R0I;
      break;
    };



  if(IPOT<4) _RCCC[IPOT]=RC/AT3;


  RMCH=_RMCHD[I];

  if(RMCH==0)   RMCH=max(RR+7.*AR, RI+7.*AI);

  NPINT=_NPD[I];
  NP=NPINT-1;
  ANP=NP;
  DR=RMCH/ANP;

  ZZO=1.4400*AZP*AZT;
  ZZI=ZZO/(2.0*RC);
  //
  R=0.0;


  ///========================= start cikl L19
  for(I=1; I<=NPINT; I++) {
      EXR=(R-RR)/AR;
      if(EXR> 40) VV0=0;
      else if(EXR<-40) VV0=-V0;
      else             VV0=-V0/(1.0+exp(EXR));

      if(ZZO==0.) VC=0.;
      else if(R>RC)    VC=ZZO/R;
      else             VC=ZZI*(3.0-R*R)/(RC*RC);

      VR[I]=VV0+VC;
      //---------------
      if(IMG!=1)    {
          EXI=(R-RI)/AI;

          if(EXI> 40 ) VW0=0;
          else if(EXI<-40 ) VW0=-W0;
          else         VW0=-W0/(1.0+exp(EXI));

          VI[I]=VW0;
        }
      else {
          EXI=(R-RI)/AI;
          if(EXI>40 || EXI<-40) VI[I]=0.;
          else {
              EXI=exp(EXI);
              VI[I]=-4.0*W0*EXI/pow((1+EXI),2);
            }
        }
      //---------------

      R+=DR;
    }
  ///========================= end cikl L19




  AM=AMP*AMT/(AMP+AMT);
  //
L10001:
  AK=0.218898*sqrt(AM*ECM);
  PLSQ=10.0*PI/pow(AK,2);
  EP=ECM*(AMP+AMT)/AMT;
  ETAC=0.157404*AZP*AZT*sqrt(AMP/EP);
  //L58:
  Q=AK*RMCH;
  DX=AK*DR;
  //
  LL = min(Q-ETAC+6.0, 1.34*(Q-2.0*ETAC+20.0));
  LL = min(1.5*Q+1.0, LL);
  LL = max(LL,0);

  if(LL>140) LL=140;

  ALX=1.25*(Q-ETAC)+1.0;
  LLX=ALX;
  LX=min(max(LL,LLX),142);
  //
  //L63:
  //
  //
  //***  SPHERICAL BESSEL FUNCTIONS  *******************************


  // =============================
  if(AZP==0) {

      GC[1]=-cos(Q)/Q;
      GC[2]=(GC[1]-sin(Q))/Q;

      if(Q-0.01 <=0){
          FC0=1.-Q*Q/6.;
          if(LX <= 1) {
              FC[1]=FC0;
              FC[2]=Q*(0.3333-Q*Q/30.);
              goto L3200;
            }
        }
      else   FC0=sin(Q)/Q;

      L2=LX+5;
      AL2=L2;
      FC[L2+1]=1.0E-10;
      FC[L2]=1.0E-10*(2.0*AL2+1.0)/Q;
      L3=L2-1;

      for(LL1=1; LL1<=L3; LL1++) { //L3106
          L1=L2-LL1;
          FL1=L1;
          FC[L1]=(2.0*FL1+1.0)*FC[L1+1]/Q-FC[L1+2];
          if(FC[L1]-1.0E+30 > 0)
            {
              for(L4=L1; L4<=L2; L4++) FC[L4] = 1.0e-10 * GC[L4];
              break;
            }
        }


      ZZ=FC0/FC[1];
      for(L1=1; L1<=2; L1++) FC[L1]*=ZZ;

      if(LX-1 > 0 ) {
          L2 -= 4;
          for(L1=3; L1<=L2; L1++) {  //L3109
              AL1=L1;
              FC[L1]*=ZZ;
              GC[L1]=(2.0*AL1-3.0)*GC[L1-1]/Q-GC[L1-2];
            }
        }
    } //L3200 -3400 =======================================================
  else {

      //***  COULOMB WAVE FUNCTIONS  *****************************
      ETA=ETAC;
      ETASQ=ETA*ETA;
      ETA2=2.0*ETA;

      if(ETA < 1 ) RX=2.0;
      else
        {
          RX=ETA2;
          if(ETA > 2)
            {
              ETA6=pow(ETA,0.166667);
              A1=(0.495957017E-01+0.245519918E-02/ETASQ+0.253468412E-03/pow(ETASQ,2))/pow(ETA,1.333333);
              A2=(0.172826037E+00+0.358121485E-02/ETASQ+0.907396643E-03/pow(ETASQ,2))/pow(ETA,0.666667);
              B1=1.0-(0.888888890E-02+0.910895810E-03/ETASQ)/ETASQ;
              B2=1.0+(0.317460317E-03+0.311782468E-03/ETASQ)/ETASQ;
              FC[1]=0.70633264*ETA6*(B1-A1);
              GC[1]=1.22340402*ETA6*(B1+A1);
              FP= 0.408695732E+00/ETA6*(B2+A2);
              GP=-0.707881773E+00/ETA6*(B2-A2);
              goto L3403;
            }
        }

      //***  ETA.LE.2.0  *************


      SMAX=(22.0+ETA2)/RX;
      COH=0.01*SMAX;
      COA=PI*ETA2;
      COC=RX*sqrt((1.0-exp(-COA))/COA);

      COAC[6]=COAC[1]=0.3298600;
      COAC[5]=COAC[2]=1.302083;
      COAC[4]=COAC[3]=0.8680556;

      for(N=1; N<=6; N++) COE[N]=0.;

      COAS=COH=0.10;

      //--------------------------
      do {
          COAA[5]=0.;
          COAA[6]=0.;
          COAS-=COH;

          for( N=1; N<=6; N++) {
              COAB=exp(-RX*COAS+2.0*ETA*atan(COAS));
              COAY[5]=COAB;
              COAY[6]=-COAB*COAS;

              COAS+=COH;
              for(M=5; M<=6; M++) COAA[M]+=COAC[N]*COAY[M];
            }
          for(M=5; M<=6; M++) COE[M] +=COH*COAA[M];

        } while(COAS-SMAX <=0 );
      //--------------------------

      COAS=COAT=0.05;

      //--------------------------
      do {
          for(N=1; N<=4; N++) COAA[N]=0.;

          COAS-=COAT;

          for( N=1; N<=6; N++) {
              COAE=tanh(COAS);
              COAA1=COAE*RX-2.0*ETA*COAS;
              COAB=sin(COAA1);
              COAA1=cos(COAA1);
              COAD=1.0-COAE*COAE;
              COAY[1]=COAD*COAA1;
              COAY[3]=-COAD*COAB;
              COAY[2]= COAE*COAY[3];
              COAY[4]=-COAE*COAY[1];

              COAS+=COAT;
              for(M=1; M<=4; M++)  COAA[M]+=COAC[N]*COAY[M];
            }
          for(M=1; M<=4; M++) COE[M]+=COAT*COAA[M];

        } while(COAS <= 10.0);
      //--------------------------


      FC[1]=COC*COE[1];
      GC[1]=COC*(COE[3]+COE[5]);
      FP=COC*COE[2]+FC[1]/RX;
      GP=COC*(COE[4]+COE[6])+GC[1]/RX;


L3403:
      D=Q-RX;
      N=fabs(D);
      STEP = min(fabs(D),1.0);
      if(D==0  ) goto L3470;

      if(D<0) {
          STEP=-STEP;
          STP=0.5;
        }
      else    STP=1.0;


      if (N>0) {
          IFN=0;
          MM=1;

L3450:        if(RX-2.0 <=0) {
              STEP*=STP;
              STP=1.0;
            }
        }

L3456:
      SIG[1]=GC[1];
      SIG[2]=STEP*GP;
      SIG[3]=STEP*STEP*(ETA2-RX)*SIG[1]/2.0/RX;
      GC[1]=SIG[1]+SIG[2]+SIG[3];
      GP=SIG[2]+2.0*SIG[3];

      for( J=4; J<=50; J++) { //L3458
          FQ=J-2;

          SIG[J]=(
                   STEP*(FQ*FQ-FQ)*SIG[J-1]+
                 STEP*STEP*(RX-ETA2)*SIG[J-2]+
                 STEP*STEP*STEP*SIG[J-3]
                 )
              /
              (-RX*FQ*(FQ+1.0));


          GC[1]+=SIG[J];
          GP+=(FQ+1.0)*SIG[J];
          if(fabs(SIG[J]/GC[1]) <= 1.0E-8 &&
             fabs(SIG[J]*FQ/GP) <= 1.0E-8) break;
        } //continue;L3458:


      GP/=STEP;
      RX+=STEP;

      if (IFN <= 0) {

          if (IFN == 0) {
              MM++;
              if(MM<=N)goto L3450;
            }

          IFN=1;

          if(fabs((RX-Q)/Q)>1.0E-6) {

              STEP=Q-RX;

              if (STEP==0) goto L3470;
              if (STEP < 0 )
                if(RX-2.0 <= 0) {
                    if(STEP<-.5)STEP=-.5;
                    IFN=-1;
                  }
              goto L3456;
            }
        }

L3470: GC[2]=((1.0/Q +ETA)*GC[1]-GP)/sqrt(1.0+ETASQ);
      L1 = LX+1;

      if (L1>2)
        for(I=3; I<=L1; I++) {
            FL=I-2;
            GC[I]=((2.0*FL+1.0)*(ETA+FL*(FL+1.0)/Q)*GC[I-1]-(FL+1.0)*sqrt(pow(FL,2)+ETASQ)*GC[I-2])/FL/sqrt(pow((FL+1.0),2)+ETASQ);
          }

      L5=L1+5;
      FC[L5+2]=0.0;
      FC[L5+1]=1.0E-10;


      for(I=1; I<=L5; I++) { //L3478
          I1 = L5-I;
          FL = I1+1;
          FC[I1+1]=((2.0*FL+1.0)*(ETA+FL*(FL+1.0)/Q)*FC[I1+2]-FL*sqrt(pow((FL+1.0),2)+ETASQ)*FC[I1+3])/(FL+1.0)/sqrt(pow(FL,2)+ETASQ);
          if (FC[I1+1]>=1.E14)
            for(I2=I1; I2<=L5; I2++) FC[I2+1]=1.0E-22*FC[I2+1];
        } //continue; L3478:


      RENORM=sqrt(1.0+ETASQ) *(FC[1]*GC[2]-FC[2]*GC[1]);

      for(I1=1; I1<=L5; I1++) FC[I1]/=RENORM;

      FP=(1.0+FC[1]*GP)/GC[1];

      // L3485:
      //      ETAX = ETASQ +16.0;

      //      SIG0=-ETA+ETA/2.0*log(ETAX)+3.5*atan(ETA/4.0)-(atan(ETA)+atan(ETA/2.0)+atan(ETA/3.0))-ETA/12.0/ETAX*(1.0+(ETASQ-48.0)/30.0/pow(ETAX,2)+(pow(ETA,4)-160.0*ETASQ+1280.0)/105.0/pow(ETAX,4));

    }  ////  ---------------------- 3400-3200

L3200:

  //
  SIGR=0.0;
  Li=0;

L23:  FL=Li;

  if(AZP!=0.) {
      FSUBL=FC[Li+1];
      GSUBL=GC[Li+1];
      if(Li==0) {
          GLP=GP;
          FLP=FP;
        }
      else   {
          TEMPA=sqrt(FL*FL+pow(ETAC,2));
          TEMPB=FL*FL/Q+ETAC;
          GLP=(TEMPA*GC[Li]-TEMPB*GC[Li+1])/FL;
          FLP=(TEMPA*FC[Li]-TEMPB*FC[Li+1])/FL;
        }
    }
  else {
      FSUBL=Q*FC[Li+1];
      GSUBL=-Q*GC[Li+1];
      if(Li==0)
        {
          GLP=-FSUBL;
          FLP=GSUBL;
        }
      else
        {
          GLP=-Q*GC[Li]+FL*GC[Li+1];
          FLP=Q*FC[Li]-FL*FC[Li+1];
        }
    }

  SSUBL=Q/(FSUBL*FSUBL+GSUBL*GSUBL);
  DELSL=(FSUBL*FLP+GSUBL*GLP)*SSUBL;
  //
  //***  INTEGRATION OF SCHROEDINGER EQUATION IN WELL  *****************
  //
  for(I=1; I<=500; I++) {
      URE[I]=0;
      UIM[I]=0;
    }

  RXX = FL*pow(0.1,(50.0/max(FL,1.0)));
  II = RXX/DX+1.0;

  if (II-1 > 0) {
      URE[II] = 1.0;
      XII=II;
      XII=XII/(XII-1.0);
      URE[II+1]=pow(XII,FL);
      UIM[II] = URE[II];
      UIM[II+1]=URE[II+1];
    }
  else  {
      II=1;
      URE[II]=0.;
      URE[II+1]=1.0E-8;
      UIM[II]=0.;
      UIM[II+1]=1.0E-8;
    }

  XX=II;
  XX*=DX;
  I=II;
  DX2=DX*DX;
  //4121 DRE=2.0-DX**2+FL*(FL+1.0)*DX**2/XX**2+VR(I+1)*DX**2/ECM

  do {
      DRE=2.0+(-1. +FL*(FL+1.0)/(XX*XX)+VR[I+1]/ECM)*DX2;
      DI=VI[I+1]*DX2/ECM;
      URE[I+2]=DRE*URE[I+1]- DI*UIM[I+1]-URE[I];
      UIM[I+2]=DRE*UIM[I+1]+ DI*URE[I+1]-UIM[I];

      if(fabs(URE[I])>1.E16) {
          MI=I+1;
          for(M=1; M<=MI; M++) {
              I=M-1;
              URE[I+2]*=1.0E-22;
              UIM[I+2]*=1.0E-22;
            }
        }

      XX+=DX;
      I++;

    } while (I<=NP);


  //
  //***  NUCLEAR PHASE SHIFTS AND TRANSMISSION COEFFICIENTS  **************
  //
  K1=NPINT+1;
  L1=Li+1;
  if(max(fabs(URE[K1-1]),fabs(UIM[K1-1]))-1.0E16) {;}//3705,3705,3701; //MPK What was this supposed to do???
  //L3701:
  for(I=1; I<=K1; I++)
    {
      URE[I]=1.0E-22*URE[I];
      UIM[I]=1.0E-22*UIM[I];
    }
  //L3705:
  DDD=pow(URE[K1-1],2)+pow(UIM[K1-1],2);
  DEN=2.0*DDD/ANP;
  FLR =(URE[K1-1]*(URE[K1]-URE[K1-2])+UIM[K1-1]*(UIM[K1]-UIM[K1-2])) /DEN;
  FLI=(-UIM[K1-1]*(URE[K1]-URE[K1-2])+URE[K1-1]*(UIM[K1]-UIM[K1-2])) /DEN;
  DEN = (pow(GSUBL,2)+pow(FSUBL,2))*(pow((FLR-DELSL),2)+pow((FLI-SSUBL),2));
  ETARE  =   ((pow((FLR-DELSL),2)+pow(FLI,2)-pow(SSUBL,2))*(pow(GSUBL,2)-pow(FSUBL,2)) +4.0*FSUBL*GSUBL*(FLR-DELSL)*SSUBL)/DEN;
  ETAIM  =   2.0*((pow(GSUBL,2)-pow(FSUBL,2))*SSUBL*(FLR-DELSL)-FSUBL*GSUBL* (pow((FLR-DELSL),2)+(pow(FLI,2)-pow(SSUBL,2))))/DEN;
  TL[L1] = 1.0-(pow(ETARE,2)+pow(ETAIM,2));
  SIGL[L1]=PLSQ*(2.0*FL+1.0)*TL[L1];
  SIGR+=SIGL[L1];

  Li++;
  if(Li<=LL)goto L23;

}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//
//   for Quantum-Mechanical
//   get T0
//
double GetTCN0(double Ecm, double V0, double omega)
{
  double t =  - 2.* PI / omega * (Ecm-V0);
  double TCN0 = 0;

  if(t < 100)
    TCN0 = 1. /  (1.+exp(t));

  return TCN0;
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//   for Quantum-Mechanical
//
//   get Energy from M
//
double GetEnergyFromL(double Ecm, int L, double reduced_mass, double Rf)
{
  double t =  2.* aem_MeV * reduced_mass * Rf * Rf;
  if(t<=0) return 0;

  double E = double(L*(L+1))*hc_MeV_fermi*hc_MeV_fermi/ t;

  return Ecm-E;
}

