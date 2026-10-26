#include <QByteArray>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "pace_random.h"

static void require(bool ok, const char* message) {
    if(!ok) throw std::runtime_error(message);
}
static void seed(const char* value) {
    qputenv("PACE4_SEED", value);
    require(InitRandom(), "Initialize fixed seed");
}
int main() {
    try {
        seed("13");
        std::vector<double> expected;
        for(int i=0;i<1000;++i) {
            const double value=MyRandom();
            require(value>0. && value<1., "Random values must be strictly inside (0,1)");
            expected.push_back(value);
        }
        const long legacyNumerators[] = {521016965,1441301936,346100192,1520210868,1571110117};
        for(int i=0;i<5;++i)
            require(std::abs(expected[i]*2147483647.-legacyNumerators[i])<1.e-6,"Legacy recurrence compatibility");
        seed("13");
        for(double value:expected) require(MyRandom()==value, "Explicit seed must replay exactly");
        require(CurrentRandomSeed()==13, "Report external seed exactly");
        const double equal[] = {0.,1./3.,2./3.,1.};
        seed("13");
        const int selected=SampleInitialSpin(equal,3);
        require(selected==(expected[0]<=1./3.?1:expected[0]<=2./3.?2:3), "CDF threshold lookup");
        require(MyRandom()==expected[1], "Spin selection must consume exactly one draw");
        int counts[4]={};
        seed("24681357");
        constexpr int samples=300000;
        for(int i=0;i<samples;++i) ++counts[SampleInitialSpin(equal,3)];
        for(int i=1;i<=3;++i) require(std::abs(double(counts[i])/samples-1./3.)<0.005,"Equal weights must produce equal bins");
        std::cout << "PASS one-draw CDF sampling: " << counts[1] << ' ' << counts[2] << ' ' << counts[3] << '\n';
        const double sparse[] = {0.,0.,0.,1.};
        for(int i=0;i<1000;++i) require(SampleInitialSpin(sparse,3)==3,"Zero-probability bins must remain empty");
        for(const char* special:{"123459876","2024023771","1","2147483646"}) {
            seed(special);
            for(int i=0;i<1000;++i) { const double value=MyRandom(); require(value>0.&&value<1.,"Degenerate seed must not lock generator at zero"); }
        }
        for(const char* invalid:{"0","-1","2147483647","nonsense","1.5"}) {
            qputenv("PACE4_SEED",invalid);
            require(!InitRandom(),"Invalid seed must be rejected");
        }
        qunsetenv("PACE4_SEED");
        require(InitRandom(),"Entropy seed initialization");
        require(CurrentRandomSeed()>0&&CurrentRandomSeed()<2147483647u,"Entropy seed range");
        std::cout << "PASS exact replay, seed validation, nonzero states, sparse distribution\n";
    } catch(const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n'; return 1;
    }
    return 0;
}
