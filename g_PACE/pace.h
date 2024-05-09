#ifndef PACE_H
#define PACE_H

#include <QTextStream>
#include "ParticleStates.h"

//WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

class PACE
{
public:
    PACE(){};

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
            int *MAXJS, double *EBIN, int IXMIN, QTextStream& s);
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
