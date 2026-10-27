#include <QApplication>
#include <QTextStream>
#include <QString>
#include <QTemporaryDir>
#include "ParticleStates.h"
#define private public
#include "pace.h"
#undef private
#include "pace_version.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>
#include <vector>

QString FileArg;
int fontsizeGlobal = 9;
int useHighDpiScaling = 1;
bool PACESilentMode = true;
extern FILE* f09;
extern double _FYRST, _FACLA, _EROT[], _DISCPR, _PROB[], _BE[], _FBARR[], _BARMAX;
extern double _FGE1, _FGM1, _FGE2, _FGM2, _EMIN[], _EMAX[], _VQ0[], _VQS[];
extern int _MAXC, _IPROB, _MJ[], _cl_K, _IDIST, _IZR[], _INR[];
extern void TLLL(double, int, int, int, double[25], int&, QTextStream&);
extern void TLOM(int, double, double, double, double, double, int&, double[Max_MOM+1]);

extern double _EX, _TT, _ERGMIN, _BARFAC;
extern int _IADEF, _NDISC, _NOSHL, _IS;
extern void GCLVD(double&, double, double, double, double, double, double, double,
                  int&, int&, int, int*, int*, double*, double*, QTextStream&);
extern void FIND(double, int&, double*, int, int, QTextStream&);

static void require(bool ok, const char* message) {
    if(!ok) throw std::runtime_error(message);
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    PACE pace;
    QString output;
    QTextStream stream(&output);
    QTemporaryDir temporary;
    require(temporary.isValid(), "Create temporary diagnostics directory");
    f09 = std::fopen(temporary.filePath("physics.rtf").toLocal8Bit().constData(), "w");
    require(f09 != nullptr, "Open temporary diagnostics file");
    int failures = 0;
    auto test = [&](const char* name, auto body) {
        try { body(); std::cout << "PASS " << name << '\n'; }
        catch(const std::exception& e) { ++failures; std::cerr << "FAIL " << name << ": " << e.what() << '\n'; }
    };
    test("all nonzero-FYRST rotational levels and repeated calls", [&]() {
        _MAXC = 21;
        _FACLA = 10;
        for(double factor : {1., .625, 0.}) {
            _FYRST = factor;
            std::fill(_EROT, _EROT + Max_MOM + 1, -12345.);
            double f=0,a=0,sigma=0,fact=0,delta=0,aconst=0,sr=.5;
            pace.GCLDP(59, 28, 0, f, a, sigma, fact, delta, aconst, sr, 1, stream);
            for(int j=1; j<=_MAXC+9; ++j) {
                const double spin=j-1+sr;
                const double expected=factor ? pace.YRAST(59,28,spin)*factor : (spin+.5)*(spin+.5)/(2*sigma);
                require(std::abs(_EROT[j]-expected) <= 1e-12*std::max(1.,std::abs(expected)), "Rotational level missing or stale");
            }
        }
    });
    test("half-integer dipole transition and forbidden zero-to-zero", [&]() {
        int maxj[3]={0,1,0}, maxjs[3]={0,1,0};
        double ebin[3]={0,9,0}, rlev[3]={0,1,0};
        _FGE1=1; _FGM1=_FGE2=_FGM2=0;
        _DISCPR=0; _BARMAX=30; _FBARR[1]=100;
        fusion_event c; c.init(28,31,10,1,0);
        double probm=0;
        _IPROB=0;
        pace.CHNPRB(4,1,c,.5,0,rlev,probm,maxj,maxjs,ebin,1,1,stream);
        require(_IPROB==1 && _PROB[1]>0, "J=1/2 to J=1/2 E1 transition was discarded");
        const double expected=6.28318*std::pow(59.,.666667)*6.8e-8;
        require(std::abs(_PROB[1]/expected-1.)<1e-12, "Incorrect dipole transition strength");
        c.init(28,32,10,1,0); probm=0; _IPROB=0;
        pace.CHNPRB(4,1,c,0,0,rlev,probm,maxj,maxjs,ebin,1,1,stream);
        require(_IPROB==0, "J=0 to J=0 E1 transition must remain forbidden");
    });
    test("optical transmission coefficients against original Fortran", [&]() {
        // Original PACE3.for TLOM and OPTPOT, unchanged, compiled with
        // GNU Fortran11.2 -fdefault-real-8 -fdefault-double-8 to match double.
        // Initial optical defaults use OPTPOT(channel,28,32), as in PACE3.
        const double expected[3][6]={
            {.77123367599980031,.69420494736765226,.74079179488658464,
             .87825484834474044,.46506178923713037,.54258313572146100},
            {.99449826743775360,.64005841032247224,.92275049760470063,
             .38866673007472952,.20579068887490148,.030455763122249624},
            {.99739894598404377,.99775788738702742,.99709329251273782,
             .99481151781558597,.99569012400522350,.98301883111167632}
        };
        const int expectedL[3]={10,10,18};
        for(int channel=1; channel<=3; ++channel) {
            pace.OPTPOT(channel,28,32);
            double tl[Max_MOM+1]={}; int ll=0;
            TLOM(channel,channel==3?4.:1.,channel-1.,channel==3?56.:59.,
                 channel==3?26.:27.,channel==3?15.:10.,ll,tl);
            require(ll==expectedL[channel-1], "Optical partial-wave limit differs from Fortran");
            double maxError=0;
            for(int l=0; l<6; ++l)
                maxError=std::max(maxError,std::abs(tl[l+1]-expected[channel-1][l]));
            std::cout << std::setprecision(17) << "Optical channel " << channel
                      << " maximum absolute transmission difference: " << maxError << '\n';
            // Agreement to six decimal places allows numerical-integration and
            // compiler differences while detecting the wrong interior potential.
            require(maxError<1e-6, "Optical transmission differs from Fortran reference");
        }
    });
    test("charged channel after an empty residual-energy bin", [&]() {
        _IDIST=0; _cl_K=0;
        for(int channel=1; channel<=3; ++channel) {
            _EMIN[channel]=.01; _EMAX[channel]=40;
            pace.OPTPOT(channel,28,32);
        }
        double tl[25]={}; int ll=0;
        TLLL(0,0,28,32,tl,ll,stream);
        for(int channel=1; channel<=3; ++channel) _VQS[channel]=_VQ0[channel];
        _IZR[2]=27; _INR[2]=32; _BE[2]=0;
        TLLL(10,2,27,32,tl,ll,stream);
        require(tl[1]>=1e-10, "Test requires an open proton optical channel");
        int maxj[3]={0,1,1}, maxjs[3]={0,1,2};
        double ebin[3]={0,5,10}, rlev[3]={0,1,0};
        fusion_event c; c.init(28,32,20,1,0);
        double probm=0; _DISCPR=0; _IPROB=0;
        pace.CHNPRB(2,2,c,0,.5,rlev,probm,maxj,maxjs,ebin,1,1,stream);
        require(_IPROB==1 && _PROB[1]>0, "Search stopped before the lower, populated residual-energy bin");
        require(_MJ[1]==79696+1, "Wrong lower-bin transition selected");
    });
    test("fission barriers cover every supported initial spin", [&]() {
        _MAXC=Max_MOM; _FYRST=0; _FACLA=10; _NOSHL=1; _NDISC=0; _IADEF=1;
        _IS=0; _BARFAC=1; _BARMAX=30; _ERGMIN=143.076;
        _EMIN[4]=0; _EMAX[4]=20;
        std::fill(_FBARR, _FBARR+Max_MOM+1, -12345.);
        std::vector<int> maxj(DIM_EBIN), maxjs(DIM_EBIN);
        std::vector<double> ebin(DIM_EBIN), rlev(DIM_RLEV);
        double masses[5]={}; int mebin=0, ixpr=0;
        fusion_event c; c.init(70,102,143.076,149,0);
        pace.AMASS(masses,c,4,mebin,maxj.data(),maxjs.data(),ebin.data(),rlev.data(),ixpr,_BE,stream);
        for(int j=1; j<=Max_MOM; ++j) {
            const double expected=pace.SIERK(172,70,j-1,1);
            require(std::isfinite(_FBARR[j]) && std::abs(_FBARR[j]-expected)<1e-12,
                    "Supported high-spin barrier was not initialized");
            require(_FBARR[j]>=0, "Nonnegative tabulated barrier became negative");
        }
    });
    test("level densities include all supported spins", [&]() {
        _MAXC=Max_MOM; _IADEF=2; _EX=0; _TT=1;
        std::fill(_EROT,_EROT+Max_MOM+1,0.);
        std::vector<int> maxj(DIM_EBIN), maxjs(DIM_EBIN);
        std::vector<double> ebin(DIM_EBIN), rlev(DIM_RLEV+Max_MOM, -12345.);
        double alit=20; int kk=0, mebin=0;
        GCLVD(alit,1,1,0,1.5,2,1,0,kk,mebin,0,maxj.data(),maxjs.data(),ebin.data(),rlev.data(),stream);
        require(mebin==1 && maxj[1]==Max_MOM && kk==Max_MOM,
                "Level-density spin range still stops at the legacy 110-spin limit");
        require(rlev[Max_MOM]>0 && std::isfinite(rlev[Max_MOM]), "Highest supported spin has no density");
    });
    test("level-density capacity is checked before writing", [&]() {
        _MAXC=101; _IADEF=2; _EX=0; _TT=1;
        std::fill(_EROT,_EROT+Max_MOM+1,0.);
        std::vector<int> maxj(DIM_EBIN), maxjs(DIM_EBIN);
        std::vector<double> ebin(DIM_EBIN), rlev(DIM_RLEV+Max_MOM,-12345.);
        int kk=DIM_RLEV-100, mebin=0; double alit=20;
        GCLVD(alit,1,1,0,1.5,2,1,0,kk,mebin,0,maxj.data(),maxjs.data(),ebin.data(),rlev.data(),stream);
        require(alit==0 && mebin==0, "Insufficient level storage was not rejected before committing a bin");
        require(std::all_of(rlev.begin()+DIM_RLEV-99,rlev.end(),[](double v){ return v==-12345.; }),
                "Level-density capacity check occurred after writing the overflowing row");
    });
    test("FIND respects the last valid energy bin", [&]() {
        _IZR[1]=70; _INR[1]=101;
        double ebin[]={0,0,1,2,1.5,4}; int index=0;
        FIND(2,index,ebin,3,1,stream);
        require(index==3, "FIND scanned beyond the final valid bin");
        FIND(1.75,index,ebin,3,1,stream); require(index==2,"Wrong interior energy bin");
        FIND(0,index,ebin,3,1,stream); require(index==1,"Wrong first energy bin");
    });
    test("fission interpolation at the final tabulated rotational row", [&]() {
        for(double mass : {60.,172.})
            for(int spin=150; spin<=199; ++spin)
                require(pace.FISROT(mass,mass==60?28:70,spin,1)==0,
                        "Barrier at the zero-height table edge must remain zero");
    });
    std::fclose(f09); f09=nullptr;
    return failures ? 1 : 0;
}
