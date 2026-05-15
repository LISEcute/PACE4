#include <iostream>
#include <QDebug>
#include <string>
#include <QDataStream> // for binary i/o
#include <QFile>
#include <QIODevice>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))
#include "ftype.h"
#include "ParticleStates.h"

#define TWO_WORDS(Hw, Lw) ((int)(((unsigned int)(Lw&0xFFFF)) + (Hw<<16)))

extern FILE *f09;
extern FILE *f_particles;

char s10m2[5] = "    ";
extern char *s10star;
extern char *s40star;
extern double mzsqrt(double X);
void STATIS();
void QQQQ(int K, int IB, QTextStream &s);

int  NDIST[31][ 21][7], IEBIN[31],JBIN[21];
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
void STATIS(fusion_event **a, const char* filename_evt,
            ClassParticleFlags &ParticleFlags, QTextStream& s)
{
  //     WRITTEN BY A. GAVRON

  double  ECM[6], DSPIN[6],  PF[31];
  //int    NSPC[31][3], NPF[31], NTFISS[31][21][7];
  //int    NSPC[31][4], NPF[31], NTFISS[32][22][8];
  int    NSPC[32][4], NPF[31], NTFISS[32][22][8];
  //=========      common/QMMJ/
  extern double /*_FYRST,_FACLA,*/_EROT[Max_MOM+1];
  //extern int    _MAXC;


  //-------    common /XQSIG/
  //extern double  _SIGMA;
  extern int _NCASC;

  //=========      common/XQRES/
  extern double _ENERGY;//,_VZC;
  extern int   _ITRAC;
  int IFISS;
  int Nnpa[4]= {0,0,0,0};
  //=========      common/OUT/

  //extern double _ERGCS[5017];
  //extern int    _IZCS[5017],_INCS[5017],_JCS[5017],_MJCS[5017];

  int IZQ[201], INQ[201],  NTOT[6];
  double YJM[31][21];

  const char *NTITL[7]={"    ","NEUT","PROT","ALPH","G-E1","G-E2","FISS"};
  //QString NTITL[7]={"    ","NEUT","PROT","ALPH","G-E1","G-E2","FISS"};
  const char *particle_format ="%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%-7.1e\t%-6.1f\t%-6.1f\t%-6.1f\t ";
  const char *particle_title1 = "Decay\tN\tN\tChain\tZ_f\tN_f\tZc\tN_c\tJ_c\tJf\tM_Jc\tFission\tEx_i\tEx_f\tEp_Lab\tAp_Lab\n";
  const char *particle_title2 = "Mode\tmode\tAll\t\tfinal\tfinal\temitter\temitter\tinit\tfinal\tproj-n\tprob\tMeV\tMeV\tMeV\tdeg\n";

  /*              IBUF2[0]=HIGH_WORD(I);
              IBUF2[1]=LOW_WORD(I);
              IBUF2[2]=(short int)(MODE+IE12);
              IBUF2[3]=(short int)csi.J;
              IBUF2[4]=(short int)csf.J;
              IBUF2[5]=(short int)csi.MJ;
              IBUF2[6]=(short int)(c.Ex+1.);
              IBUF2[7]=(short int)(EP*10.+1.);
              IBUF2[8]=(short int)c.Z;
              IBUF2[9]=(short int)c.N;

      MODE = IBUF2[2];
      JC   = IBUF2[3];
      JF   = IBUF2[4];
      MP   = IBUF2[5];
      IE   = IBUF2[6];
      IP   = IBUF2[7];
      IZC  = IBUF2[8];
      INC  = IBUF2[9];
*/

  //      THE ABOVE CAN PROVIDE DETAILED INFORMATION ON ALL CASCADES
  //      LEADING TO ANY SPECIFIC FINAL NUCLEUS.
  //      INDX -  IZCS(INDX) AND INCS(INDX) ARE THE Z AND N OF THE FINAL
  //       NUCLEUS IN THE CHAIN OF DECAY.
  //      MODE - TYPE OF PARTICLE EMITTED. 1-N, 2-P, 3-ALPHA, 4-GAMMA
  //      JC, JF - INITIAL AND FINAL SPIN INDICES FOR THIS PARTICLE
  //       EMISSION. (EACH READ STATEMENT CORRESPONDS TO 1 EMITTED PART.)
  //      MP - THE PROJECTION OF JC ON THE Z AXIS. ( FRACTIONAL SPINS
  //       ARE NEGLECTED FOR THE PROJECTION. FOR JC,JF, THE ACTUAL
  //       SPIN = JC,JF-1   IN EVEN MASS NUCLEUS
  //            = JC,JF -1/2   IN ODD MASS NUCLEUS )
  //      IE - EXCITATION ENERGY AT EMITTING LEVEL
  //      IP - PARTICLE ENERGY (C.M.) MULTIPLIED BY FACTOR OF  2(OLD)  -- really 10 now
  //      IZC, INC - Z N , OF EMITTING NUCLEUS
  //      IF IZC IS NEGATIVE,IT MEANS THAT THE NUCLEUS FISSIONED
  //

  int I,K,L,J,NZN, IZPR;
  double SFE,TFE,SFS,TFS,GFS,SFM,TFM, Fprob;
  double AM,EC, DSP, SPIN;
  double SE,TE,SS,TS,EXEN;
  double Ex_i, Ex_f, a_Angle, a_Elab;

  int INDX, MODE, JC, JF, MP, IE, IP, JX, IB, INUC, IE1, JC1,IEE1, NI, IZC, INC;//, MJ;
  //int while_true_it = 0; // MPK
  //short int IBUF2[10] = {0};
  //float BUF3[5] = {0};

  double SUM=0, SUMPF=0;
  int NTT=0;
  //s10m="           ";

  char ss520[200]= "\\par  Excitation energy window - average =\\b  %-6.1f\\b0    FMHW = \\b %-6.1f\\b0 "
                   "\\par  Spin window              - average =\\b  %-6.1f\\b0    FMHW = \\b %-6.1f\\b0 ";
  //QString s520 = "\\par  Excitation energy window - average =\\b  %-6.1f\\b0    FMHW = \\b %-6.1f\\b0 \n\\par  Spin window              - average =\\b  %-6.1f\\b0    FMHW = \\b %-6.1f\\b0 ";
  char ss520b[200] ="\\par  Average fabs projection   \\b %.1f\\b0     Average rms proj \\b %.1f\\b0 ";
  //QString s520b = "\\par  Average fabs projection   \\b %.1f\\b0     Average rms proj \\b %.1f\\b0 ";

  //char ss520[200] = " Excitation energy window - average =  %-6.1f   FMHW = %-6.1f Spin window   - average = %-6.1f    FMHW =  %-6.1f ";

  //for(I=1; I<=31; I++) for(K=1; K<=3;  K++)   NSPC[I][K]=0;
  //for(I=1; I<=31; I++) for(K=1; K<=21; K++) for(L=1; L<=7; L++)   NTFISS[I][K][L]=0;

  for(I=1; I<=31; I++){
      for(K=1; K<=3;  K++){ //for(K=1; K<=3;  K++){
          NSPC[I][K]=0;
        }
    }
  for(I=1; I<=31; I++){
      for(K=1; K<=21;  K++){
          for(L=1; L<=7; L++){
              NTFISS[I][K][L]=0;
            }
        }
    }

  fprintf(f09,"\\par\\par\\par\\b\\f2\\fs25\\cf1 ***** Complete traceback diagnostic of particle and gamma emission *****\\par\\b0\\f0\\fs16\\cf0 ");
  s << "<p>&nbsp;</p><h3>***** Complete traceback diagnostic of particle and gamma emission *****</h3>";

  SFE=TFE=SFS=TFS=GFS=SFM=TFM=0;

  for(I=1; I<=30; I++) {
      PF[I]=NPF[I]=0;
      for(J=1; J<=20; J++) {
          YJM[I][J]=NDIST[I][J][6]=0;
        }
    }

  for(I=1; I<=31; I++) _EROT[I]=0.1*I-0.05;

  //      rewind(f05);

  int k=0;
  //int size;


  /*  Qt-Oleg
QString filename_evt2 = QString::fromUtf8(filename_evt) + "2";
if(QFile::exists(filename_evt2)){   QFile::remove(filename_evt2);}

if(!QFile::copy(filename_evt,filename_evt2)){qDebug()<< "copy not successfull!";}
*/

  QFile f5(filename_evt);  // Qt-Oleg
  if(!f5.open(QIODevice::ReadOnly)){
      qDebug() << "Could not open " << filename_evt << "!";
      return;
    }


  if(f_particles) {fprintf(f_particles,particle_title1); fprintf(f_particles,particle_title2);}

  //=========================

  QDataStream in(&f5);
  in.setByteOrder(QDataStream::LittleEndian);
  in.setFloatingPointPrecision(QDataStream::SinglePrecision);
  f5.seek(0);
  short int IBUF2a[10];
  //qint16 IBUF2a[10];
  float BUF3a[5];
  while(true) {

      //  size=_read(f05,IBUF2,sizeof(short int)*10);      if(size<=0)break;
      //  size=_read(f05,BUF3,sizeof(float)*5);            if(size<=0)break;

      if(!f5.bytesAvailable()){
          break;
        }
      for(int i=0;i<10;i++){
          in >> IBUF2a[i];
        }
      for(int i=0;i<5;i++){
          in >> BUF3a[i];
        }
      k++;
      INDX = TWO_WORDS(IBUF2a[0],IBUF2a[1]);

      MODE = IBUF2a[2];
      JC   = IBUF2a[3];
      JF   = IBUF2a[4];
      MP   = IBUF2a[5];
      IE   = IBUF2a[6];
      IP   = IBUF2a[7];
      IZC  = IBUF2a[8];
      INC  = IBUF2a[9];

      Fprob    = BUF3a[0];
      Ex_i     = BUF3a[1];
      Ex_f     = BUF3a[2];
      a_Elab   = BUF3a[3];
      a_Angle  = BUF3a[4];



      if(f_particles && (MODE>0 && MODE<6))
        {
          int mmode = MODE;
          if(MODE==5) mmode--;
          short int mask = 1 << (mmode-1);

          if(ParticleFlags.IsParticleStateSet(mask) &&
             ParticleFlags.CheckGate(a[INDX]->Z,a[INDX]->N))
            {
              Nnpa[3]++; Nnpa[mmode-1]++;
              fprintf(f_particles,particle_format,MODE,
                      Nnpa[mmode-1],Nnpa[3],INDX,a[INDX]->Z,a[INDX]->N,IZC,INC,
                  JC,JF,MP,
                  Fprob > 1e-10 ? Fprob : 0,Ex_i,Ex_f,a_Elab);

              if(mmode!=4) fprintf(f_particles,"%05.1f",a_Angle);
              fprintf(f_particles,"\n");
            }
        }

      if (MODE<4) continue;    //N, P, ALPHA, GAMMA

      JX=MODE-3;
      IP=min(IP,30);
      NSPC[IP][JX]++;
    }

  //=========================

  fprintf(f09,"\\par\\par\\b\\cf4 %s      Components  of  gamma  spectrum\\b0\\cf0 "
              "\\par\\b %s     Energy           E1-spec  E2-spec\\b0\\par ",s10m2,s10m2);
  s << "<p>&nbsp;</p><h4 align=center>Components of gamma spectrum </h4><table align=center cellpadding=8><tr><th>Energy</th><th>E1-spec</th><th> E2-spec</th></tr>";
  for(I=1; I<=29; I++)
    if(NSPC[I][1]+NSPC[I][2]>0) {
        fprintf(f09, "\\par %s %5.2f - %5.2f         ", s10m2,_EROT[I],_EROT[I+1]);
        s << "<tr><td>" << QString::number(_EROT[I]) << " - " << QString::number(_EROT[I+1]) << "</td>";
        if(NSPC[I][1]>0) {
            fprintf(f09, "%5d    ",NSPC[I][1]);
            s << "<td align=center>" <<  QString::number(NSPC[I][1]) << " </td>";
          } else  {           fprintf(f09, "         ");
            s << "<td></td> ";
          }
        if(NSPC[I][2]>0) {
            fprintf(f09, "%5d",NSPC[I][2]);
            s << "<td align=center>" << QString::number(NSPC[I][2]) << "</td>";
          }
        s << "</tr>";
      }

  fprintf(f09,"\\par %s  Above  %5.2f         %5d    %5d\\par",s10m2,_EROT[31],NSPC[30][1],NSPC[30][2]);
  s << "<tr><td> Above " << QString::number(_EROT[31]) << "</td><td align=center> " << QString::number(NSPC[30][1]) << "</td><td align=center>" << QString::number(NSPC[30][2]) << "</td></tr></table>";
  //printf(" %s  Above  %5.2f         %5d    %5d\n",s10m2,_EROT[31],NSPC[30][1],NSPC[30][2]);

  IB=_ENERGY*0.03333+1;
  s.flush();
  for( I=1; I<=30; I++) {
      if (I<=20) JBIN[I]=-(5*I-1);
      IEBIN[I]=I*IB;

    }

  INUC=1;

  IZQ[INUC]=a[1]->Z;
  INQ[INUC]=a[1]->N;

L60:  //lseek(f05,0,SEEK_SET);
  //qDebug() << "seek_set = " << SEEK_SET;
  //L60:  in.skipRawData(SEEK_SET);
  f5.seek(SEEK_SET);

  // Initializing the NDIST array
  for(K=1; K<=5; K++) {
      ECM[K]=DSPIN[K]=NTOT[K]=0;
      for(I=1; I<=30; I++)
        for(J=1; J<=20; J++)
          NDIST[I][J][K]=0;
    }

  if (_ITRAC==0) goto L200;


  while(true) {
      //while_true_it++;
      if(!f5.bytesAvailable()){
          break;
        }
      for(int i=0;i<10;i++){
          in >> IBUF2a[i];
        }
      for(int i=0;i<5;i++){
          in >> BUF3a[i];
        }
      k++;
      INDX = TWO_WORDS(IBUF2a[0],IBUF2a[1]);

      MODE = IBUF2a[2];
      JC   = IBUF2a[3];
      JF   = IBUF2a[4];
      MP   = IBUF2a[5];
      IE   = IBUF2a[6];
      IP   = IBUF2a[7];
      IZC  = IBUF2a[8];
      INC  = IBUF2a[9];

      //    INDX = TWO_WORDS(IBUF2[0],IBUF2[1]);

      //    MODE = IBUF2[2];
      //    JC   = IBUF2[3];
      //    JF   = IBUF2[4];
      //    MP   = IBUF2[5];
      //    IE   = IBUF2[6];
      //    IP   = IBUF2[7];
      //    IZC  = IBUF2[8];
      //    INC  = IBUF2[9];

      Fprob    = BUF3a[0];
      Ex_i     = BUF3a[1];
      Ex_f     = BUF3a[2];
      a_Elab   = BUF3a[3];
      a_Angle  = BUF3a[4];

      // size=read(f05,IBUF2,sizeof(short int)*10);    if(size<=0)break;
      // size=read(f05,BUF3,sizeof(float)*5);          if(size<=0)break;

      //    INDX = TWO_WORDS(IBUF2[0],IBUF2[1]);
      //    MODE = IBUF2[2];
      //    JC   = IBUF2[3];
      //    JF   = IBUF2[4];
      //    MP   = IBUF2[5];
      //    IE   = IBUF2[6];
      //    IP   = IBUF2[7];
      //    IZC  = IBUF2[8];
      //    INC  = IBUF2[9];

      //    Fprob    = BUF3[0];
      //    Ex_i     = BUF3[1];
      //    Ex_f     = BUF3[2];
      //    a_Elab   = BUF3[3];
      //    a_Angle  = BUF3[4];

      if (IZQ[INUC]!=a[INDX]->Z || INQ[INUC]!=a[INDX]->N) continue;

      IP=(IP+9)/10;             //     CHANGE SCALE TO 1 MEV BINS
      GFS+=Fprob;
      SFE+=IE*Fprob;
      TFE+=IE*IE*Fprob;
      SFS+=JC*Fprob;
      TFS+=JC*JC*Fprob;
      SFM+=Fprob*MP;
      TFM+=Fprob*MP*MP;
      //      IE1=min(IP/IB+1-1,30); I'm pretty sure this change will output a traceback of the neutron spectra in the COM frame, could be an option - JKS
      //      IE1=min(IE/IB+1-1,30); Proposed correction - JKS
      IE1=min(IE/IB+1,30);
      JC1=min(JC/5+1,20);
      NDIST[IE1][JC1][MODE]++;

      if (MODE<=5) {
          NTOT[MODE]++;
          ECM[MODE]+=IP;
          DSPIN[MODE]+=JC-JF;
        }
    }

  //      fprintf(f09,"\\par ");
  NZN=0;

  for(I=1; I<=_NCASC; I++)
    if (IZQ[INUC]==a[I]->Z && INQ[INUC]==a[I]->N) NZN++;

  //========================================


  for(K=1; K<=5; K++) {
      if (NTOT[K]==0) continue;

      if (IZQ[INUC]<=0) {
          IFISS=1;
          IZPR=-IZQ[INUC];
        }
      else    {
          IFISS=0;
          IZPR=IZQ[INUC];
        }

      AM=double(NTOT[K])/double(NZN);
      EC=ECM[K]/NTOT[K]-0.5;
      DSP=DSPIN[K]/NTOT[K];

      fprintf(f09,"\\par\\b  Z=%3d N=%4d  Mode=\\fc3 %.4s\\fc0    Total=%d   Out of %d events.\\b0    Multiplicity=%.2f  Average ECM=%.2f MeV"
                  " Average spin removed = %.1f",
              IZPR,INQ[INUC],NTITL[K],NTOT[K],NZN,AM,EC,DSP);
      s << "<h3> Z = " << QString::number(IZPR) << " N = " << QString::number(INQ[INUC]) << "  Mode = <span style=\"color:red\">" << QString::fromUtf8(NTITL[K]) << "</span> Total = "
        << QString::number(NTOT[K]) << "  Out of " << QString::number(NZN) << " events.  Multiplicity = " << QString::number(AM) <<
           " Average ECM = " << QString::number(EC) << " MeV Average spin removed = " << QString::number(DSP) << "</h3>";
      if (IFISS==1){
          fprintf(f09,"\\par  *************  Fission follows  particle.   Emission cascades  *************");
          s << "<p>  *************  Fission follows  particle. Emission cascades  *************</p>";
        }
      s.flush();
      QQQQ(K, IB,s);
      s.flush();
    }
  //========================================
  for(I=1; I<=_NCASC; I++) {   //L170
      for(J=1; J<=INUC; J++)
        if (IZQ[J]==a[I]->Z && INQ[J]==a[I]->N ) goto L170;

      INUC++;
      IZQ[INUC]=a[I]->Z;
      INQ[INUC]=a[I]->N;
      goto L60;
L170:
      std::cout << "L170" << std::endl;
    }



  for(K=1; K<=6; K++)
    {
      ECM[K]=DSPIN[K]=NTOT[K]=0;

      for(I=1; I<=30; I++)
        for(J=1; J<=20; J++)
          NDIST[I][J][K]=0;
    }


  //_lseek(f05,0,SEEK_SET); // mpk lseek = _lseek for ISO C++ conformant
  f5.seek(SEEK_SET);

L200:
  //size=_read(f05,IBUF2,sizeof(short int)*10);      if(size <= 0)goto L210; // mpk read = _read for ISO C++ conformant
  //size=_read(f05,BUF3,sizeof(float)*5);            if(size <= 0)goto L210; // mpk read = _read for ISO C++ conformant
  if(f5.bytesAvailable()){
      for(int i=0;i<10;i++){
          in >> IBUF2a[i];
        }
      for(int i=0;i<5;i++){
          in >> BUF3a[i];
        }
    } else {
      goto L210;
    }
  //INDX = TWO_WORDS(IBUF2[0],IBUF2[1]);
  //MODE = IBUF2[2];
  //JC   = IBUF2[3];
  //JF   = IBUF2[4];
  //MP   = IBUF2[5];
  //IE   = IBUF2[6];
  //IP   = IBUF2[7];
  //IZC  = IBUF2[8];
  //INC  = IBUF2[9];

  //Fprob    = BUF3[0];
  //Ex_i     = BUF3[1];
  //Ex_f     = BUF3[2];
  //a_Elab   = BUF3[3];
  //a_Angle  = BUF3[4];
  INDX = TWO_WORDS(IBUF2a[0],IBUF2a[1]);

  MODE = IBUF2a[2];
  JC   = IBUF2a[3];
  JF   = IBUF2a[4];
  MP   = IBUF2a[5];
  IE   = IBUF2a[6];
  IP   = IBUF2a[7];
  IZC  = IBUF2a[8];
  INC  = IBUF2a[9];



  Fprob    = BUF3a[0];
  Ex_i     = BUF3a[1];
  Ex_f     = BUF3a[2];
  a_Elab   = BUF3a[3];
  a_Angle  = BUF3a[4];

  IE1=IE/IB+1;
  IP=(IP+9)/10;
  GFS+=Fprob;
  SFE+=IE*Fprob;
  TFE+=IE*IE*Fprob;
  SFS+=JC*Fprob;
  TFS+=JC*JC*Fprob;
  SFM+=Fprob*MP;
  TFM+=Fprob*MP*MP;
  //      CHANGE SCALE TO 1 MEV BINS
  IE1=min(IE1,30);
  JC1=JC/5+1;

  JC1=min(JC1,20);

  //YJM[IE1][JC1]+=fabs(MP);
  YJM[IE1][JC1]+=abs(MP);
  if ( a[INDX]->Z < 0)
    NTFISS[IE1][JC1][MODE]++;

  PF[IE1]+=Fprob;
  NPF[IE1]++;
  if (MODE>5  || a[INDX]->Z <0) goto L200;

  NDIST[IE1][JC1][MODE]++;
  //     TRACEBACK NDIST IS RESIDUES ONLY
  NTOT[MODE]++;
  ECM[MODE]+=double(IP);
  DSPIN[MODE]+=JC-JF;
  goto L200;
L210: //fprintf(f09,"450 ");
  NZN=0;

  for(I=1; I<=_NCASC; I++)  if(a[I]->Z>=0)NZN++;

  for(I=1; I<=30; I++) {
      for(J=1; J<=20; J++) {
          SUM=NTFISS[I][J][6]+1.E-19;
          for(K=1; K<=5; K++)  SUM+=NDIST[I][J][K]+NTFISS[I][J][K];
          YJM[I][J]/=SUM;
        }
    }
  fprintf(f09,"\\par\\par\\par\\b   M states at final J vs Ex\\b0\\fs14\\par ");
  s << "<p>&nbsp;</p><h3>  M states at final J vs Ex</h3><table>";
  for(I=1; I<=30; I++) {
      IEE1=IEBIN[I]-IB;
      fprintf(f09,"\\par\\cf3  %4d - %-4d\\cf0  ",IEE1,IEBIN[I]);
      s << "<tr><td align=center style=\"blue\">" << QString::number(IEE1) << " - " << QString::number(IEBIN[I]) << "</td>";
      for(J=1; J<=17; J++){
          if(YJM[I][J] > 0) {fprintf(f09,"%6.1f",YJM[I][J]); s << "<td>" << QString::number(YJM[I][J],'f',2) << "</td>";}
          else              {
              fprintf(f09," %6.1f  \\cf2 ..\\cf0  ",YJM[I][J]);
              //s << "<td> " << QString::number(YJM[I][J]) << " .. </td>";
              s << "<td> .. </td>";
            }

        }
      s << "</tr>";
    }
  fprintf(f09,"\\par\\b    Ex / J    ");
  s << "<tr><td><b>   Ex / J</b></td> ";
  for(J=1; J<=16; J++) {
      fprintf(f09,"%6d",JBIN[J]);
      s << "<td>" << QString::number(JBIN[J],'f',2) << "</td>";
    }
  fprintf(f09,"\\fs16\\b0 ");
  s << "</tr></table>";

  const char *NTITLs[]={"    ","NEUT","PROT","ALPH","G-E1","G-E2","FISS"};
  for(K=1; K<=5; K++) {
      if (NTOT[K]==0) continue;
      AM=NTOT[K];
      AM=AM/NZN;
      EC=ECM[K]/NTOT[K]-.5;
      DSP=DSPIN[K]/NTOT[K];
      fprintf(f09,"\\par\\par\\par\\b\\fs16  Decay summary.   Mode =\\cf1  %s\\cf0\\fs14    Total = %d   Outof = %d events    "
                  "Multiplicity = %.2f\\b0 "
                  "\\par       Average Ecm = %.2f MeV   Average spin removed = %.2f\\par",NTITLs[K],NTOT[K],NZN,AM,EC,DSP);
      s << "<p>&nbsp;</p><h3> Decay summary.   Mode = " << QString::fromUtf8(NTITLs[K]) << " Total = " << QString::number(NTOT[K])
        << " Out of = " << QString::number(NZN) << " events  Multiplicity = " << QString::number(AM) << "</h3><p>"
                                                                                                        " Average Ecm = " << EC << " Average spin removed = " << DSP << "</p>";
      s.flush();
      QQQQ(K,IB,s);
    }

  fprintf(f09,"\\par\\par\\par\\fs20\\b\\cf4  Fission probability as function of excitation\\fs16\\cf0 "
              "\\par        Ex.Energy    Probability\\b0 ");
  s << "<p>&nbsp;</p><h3>Fission probability as function of excitation</h3><table><tr><th>  Ex.Energy  </th><th>  Probability</th></tr>";

  for(I=1; I<=30; I++) {
      NI=31-I;
      if (NPF[NI]<=0)continue;
      SUMPF+=PF[NI];
      PF[NI]/=NPF[NI];
      IEE1=IEBIN[NI]-IB;
      fprintf(f09,"\\par        %4d-%-4d    %10.2e",IEE1,IEBIN[NI],PF[NI]);
      s << "<tr><td align=center>" << QString::number(IEE1) << "-" << QString::number(IEBIN[NI]) << "</td><td>  " << QString::number(PF[NI],'e',2) << " </td></tr>";
    }
  s << "</table>";
  GFS = max(GFS,1e-10);
  SUMPF/=_NCASC;

  fprintf(f09,"\\par\\par\\fs16\\cf2 Total sum of fission probabilities  %.2e\\cf0 ",SUMPF);
  s << "<p>&nbsp;</p><h3 style=\"color:blue\">  Total sum of fission probabilities  " << QString::number(SUMPF, 'e',3) << "</h3>";

  SFE/=GFS;
  TFE=mzsqrt(TFE/GFS-SFE*SFE)*2.36;
  SFS/=GFS;
  TFS=mzsqrt(TFS/GFS-SFS*SFS)*2.36;
  SFM/=GFS;
  TFM=mzsqrt(TFM/GFS);

  fprintf(f09,ss520,SFE,TFE,SFS,TFS); //////////MPK CHECK
  s << "<table><tr><td>Excitation energy window - average = " << QString::number(SFE,'g',4) << "</td><td> FWHM = " <<
       QString::number(TFE,'g',3) << "</td></tr><tr><td>Spin window - average = " << QString::number(SFS,'g',3) << "</td><td> FWHM = " <<
       QString::number(TFS,'g',3) << "</td></tr></table>";
  fprintf(f09,ss520b,SFM,TFM); ////////////////MPK CHECK
  s << "<p> Average fabs projection " << QString::number(SFM,'g',3) << " Average rms proj " << QString::number(TFM,'g',3) << "</p>";
  for(I=1; I<=30; I++)
    for(J=1; J<=20; J++)
      for(K=1; K<=6; K++)
        NDIST[I][J][K]=NTFISS[I][J][K];

  for(I=1; I<=30; I++)
    for(J=1; J<=20; J++)
      NTT+=NDIST[I][J][6];

  if (NTT==0) return;

  fprintf(f09,"\\par\\b  ************************************"
              "\\par\\par\\cf3  Fission summary.\\cf0  Total events = %d "
              "\\par\\par\\cf1\\b\\fs22  Particles preceeding fission\\cf0\\fs16 "
              "\\par  ************************************\\b0 "
          , NTT);
  s << "<hr /><p>&nbsp;</p><p><b>  Fission summary.  Total events = " << QString::number(NTT)
    << "</b></p><h3><u>  Particles preceeding fission</u></h3>";

  //===================================== cikl 350 start
  for(K=1; K<=6; K++) {
      NTT=0;
      for(I=1; I<=30; I++)
        for(J=1; J<=20; J++)
          NTT+=NDIST[I][J][K];

      if (NTT==0) continue;
      fprintf(f09,"\\par\\par\\b\\fs16   Mode = \\cf3 %s\\cf0    Total number = %d\\par\\b0\\fs14 ",NTITL[K],NTT);
      s << "<p>&nbsp;</p><p> Mode = <span style=\"color:blue\">" << NTITL[K] << "  </span>Total number = " << QString::number(NTT) << "</p>";
      if(K>=6) {
          SE=TE=SS=TS=0.;

          for(I=1; I<=30; I++){
              for(J=1; J<=20; J++) {
                  SPIN=5*J-3;
                  EXEN=IEBIN[I]-0.5*IB;
                  SE+=EXEN*NDIST[I][J][6];
                  TE+=EXEN*EXEN*NDIST[I][J][6];   //check INT to double
                  SS+=SPIN*NDIST[I][J][6];
                  TS+=SPIN*SPIN*NDIST[I][J][6];
                }
            }

          SE/=NTT;
          TE=mzsqrt(TE/NTT-SE*SE)*2.36;
          SS/=NTT;
          TS=mzsqrt(TS/NTT-SS*SS)*2.36;
          fprintf(f09,ss520, SE,TE,SS,TS);
          s << "<table><tr><td>Excitation energy window - average = " << QString::number(SE,'g',4) << "</td><td> FWHM = " <<
               QString::number(TE,'g',3) << "</td></tr><tr><td>Spin window - average = " << QString::number(SS,'g',3) << "</td><td> FWHM = " <<
               QString::number(TS,'g',3) << "</td></tr></table>";
        }
      s.flush();
      QQQQ(K,IB,s);
    }
}
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

void QQQQ(int K, int IB, QTextStream& s)
{

  int ISUM, IAVR, ISGM, IEE1;
  int I,J, JSUM,L;
  double AVR,STD;
  double J5;
  s << "<table cellspacing=4>";
  for(I=1; I<=30; I++) {
      ISUM=IAVR=ISGM=0;
      for(J=1; J<=20; J++) {
          J5=5*J-3;
          IAVR+=NDIST[I][J][K]*J5;
          ISGM+=NDIST[I][J][K]*J5*J5;
          ISUM+=NDIST[I][J][K];
        }
      if (ISUM==0) continue;
      AVR=IAVR;
      AVR/=ISUM;
      STD=ISGM;
      STD=mzsqrt(STD/double(ISUM)-AVR*AVR);
      IEE1=IEBIN[I]-IB;
      fprintf(f09,"\\par %4d-%-4d  ",IEE1,IEBIN[I]);
      s << "<tr><td> " << QString::number(IEE1) << "-" << QString::number(IEBIN[I]) << "</td>";
      for(L=1; L<=16; L++){
          if(NDIST[I][L][K]>0) {
              fprintf(f09,"%6d",NDIST[I][L][K]);
              s << "<td>" << QString::number(NDIST[I][L][K]) << "</td>";
            }else  {
              fprintf(f09,"    ..");//,NDIST[I][L][K]);
              s << "<td> .. </td>";
            }
        }
      fprintf(f09,"%6d%6.1f%6.1f",ISUM,AVR,STD);
      s << "<td>"<<QString::number(ISUM)<< "</td><td>" << QString::number(AVR) <<"</td><td>" <<  QString::number(STD) << "</td></tr>";
    }
  s << "<tr><th>   Ex / J  </th>";

  fprintf(f09,"\\par\\cf2   Ex / J   ");
  for(I=1; I<=16; I++){  fprintf(f09,"%6d",JBIN[I]); s << "<td>" << QString::number(JBIN[I]) << "</td>";}
  fprintf(f09,"   sum  avrg  stdv");
  s << "<th>  sum</th><th> avrg</th><th> stdv</th></tr><tr><th> Sum </th>  ";
  fprintf(f09,"\\par\\cf3  \\b  Sum      ");
  for(J=1; J<=16; J++) {
      JSUM=0;
      for(I=1; I<=30; I++) JSUM+=NDIST[I][J][K];
      if(JSUM>0){fprintf(f09,"%6d", JSUM); s << "<td>" << QString::number(JSUM) << "</td>";}
      else    {  fprintf(f09,"    ..%6d", JSUM); s << "<td>   ..</td>";// << QString::number(JSUM);}
        }
    }
  fprintf(f09,"\\cf0\\b0");
  s << "</tr></table>";


}

