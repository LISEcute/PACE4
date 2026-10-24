#ifndef PACE_H
#define PACE_H

#include <QTextStream>
#include "ParticleStates.h"

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

class PACE
{
public:
    PACE(){}

    int PACE3(const char *filename_rtf, const char *filename_evt, const char *filename_cs,
              const char *filename_particles, const char *filename_html);
    void OPTPOT(int I,int IN, int IZ);
    ClassParticleFlags ParticleFlags;
    int _ITAKE[5];

private:

    void COMPND(int IZ, int IN, double ENERGY, int MAXC, QTextStream &s);
    void GCLDP (double AMR, double AZR, int IADEF, double &FALIT,
                double &ALIT, double &SIGSQ, double &FACT1,
                double &DELTA, double &ACONST, double &SR, double ARATIO, QTextStream& s);
    void LEVDEN(int IZ, int IN, int IFIRST, int IMAX, int *MAXJ, int *MAXJS, double *EBIN,
                int &MEBIN, double *RLEV, int &IXPR, double /*BE*/, double ARATIO, QTextStream &s);
    double YRAST(double A, double Z, double AL);
    int COMPOS(int &IZ, int &IN, double &EEXCN, double &VZC, int &MAXC,
               double &CS_CN, double &ECLOSS, int &INPUT, QTextStream& s);
    void AMASS(double *A, fusion_event &c, int MODE, int& MEBIN, int *MAXJ, int *MAXJS,
                     double *EBIN, double *RLEV, int &IXPR, double *BE, QTextStream& s);

    void CHNPRB(int IPOT,int IXR, fusion_event &c, double SC,double SE,double *RLEV,double& PROBM, int *MAXJ,
                int *MAXJS, double *EBIN, int IXPR, int IXMIN, QTextStream& s);

    double FISROT(double A, double Z, double AL, double BARFAC);
    double SIERK(double A , double Z, double AL, double BARFAC);
    void OUTEM(int ICTL, QTextStream& s, int MODE=0, double EEML=0, double AEML=0);
    void MOMENT(int ID,double A, double AP, int MODE, double EP,
                fusion_event *f, QTextStream& s);
    void TRACK(int K,int JI,int JF,double EI,double EF,double ENERGY,
               int IC,double SPROB,double ECLOSS, QTextStream& s);

};
//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#endif // PACE_H



//============================================================================

//
//
//     INPUT = 1  - PROJECTILE + TARGET
//           = 2  - COMP. NUC. , EX AND J FIXED
//           = 3  - COMP. NUC. , EX FIXED, J DISTRIBUTED
//           = 4  - COMP. NUC. ,J DISTIBUTION CALCULATED INTERNALLY
//   THE SPIN DISTRIBUTION OF THE STARTING EXCITED NUCLIDE HAS BEEN
//          DETERMINED.  THE DETERMINATION OF THE PROBABILITIES OF
//          THE PATHS OF DEEXCITATION WILL BE DETERMINED IN THE FOLLOWING
//          WAY
//         1) USING RANDOM NUMBERS, NCASC CHANNELS ARE ESTABLISHED FOR THE COMPOUND NUCLEUS.
//         2) THE HIGHEST A IN THE ARRAY IS DETERMINED.
//         3) THE FIRST NUCLEUS OF THIS A IN THE ARRAY IS DETERMINED.
//         4) THE HIGHEST ENERGY FOR THIS NUCLEUS IS DETERMINED.
//         5) THE FIRST SPIN OF THIS NUCLEUS AND ENERGY IS FOUND.
//         6) THE PROBABILITY MATRIX FOR THIS NUCLEUS ENERGY AND SPIN IS DETERMINED.
//         7) BY MEANS OF RANDOM NUMBERS, THE DEEXCITATION PATH IS CHOSEN
//            FOR ALL MEMBERS OF THE ARRAY OF THAT NUCLIDE, ENERGY AND
//            SPIN.  THE NEW Z, N, ENERGY, AND SPIN ARE SUBSTITUTED FOR
//            THE OLD VALUES IN THE ARRAY.
//         8) THE NEXT SPIN FOR THE (OLD) NUCLIDE AND ENERGY IS LOCATED.
//            IF IT EXISTS, GO TO STEP 7.
//         9) THE NEXT HIGHER ENERGY FOR THAT (OLD) NUCLIDE IS LOCATED.
//            IF IT EXISTS, GOT TO STEP 5.
//        10) THE NEXT NUCLIDE FOR THAT (OLD) A IS LOCATED.
//            IF IT EXISTS, GO TO STEP 4.
//        11) GO TO STEP 2.
//
//
//   SPROB   IS UNNORMALIZED SUM OF PROBABILITIES
//   PROB    (EQUIVALENT - PROB) IS UNNORMALIZED PROBABILITIES
//   MODES   IS NUMBER OF DIFFERENT PARTICLES (AND GAMMAS)
//           CONSIDERED FOR EVAPORATION
//   IZPART  AND  INPART  ARE Z AND N OF EMITTED PARTICLES
//     THE VECTOR ASSOCIATED WITH EACH EVENT IS AS FOLLOWS
//
//     ERGCS(J)- THE EXCITATION ENERGY
//     IZCS(J) - THE ATOMIC NUMBER ( Z )
//     INCS(J) - THE NUMBER OF NEUTRONS ( N )
//     JCS(J)  - THE SPIN INDEX  I.E.  1,2,3,   ETC.
//               SPIN = SPIN INDEX - 1  FOR EVEN A
//               SPIN = SPIN INDEX - 1/2 FOR ODD A
//     MJCS(J) - THE SPIN PROJECTION  =  0,1,2,   SPIN.
//              WE NEGLEGT 1/2 INTEGER VALUES FOR M *******
//  *** NOTE ***  THESE INTEGERS ARE   INTEGER*2  TO SAVE SPACE
//
//
//   SPIN DISTRIBUTION OF INITIAL NUCLIDE NORMALIZED TO SUM=1.
//     FOR J.GT.100 THERE IS NO COMPOUND NUCLEUS IN PRACTICE.
//     ALEV  IS PARTIAL XSECTION FOR SPIN INDEX I
//           (FROM  COMPOSE  THRU COMMON)


