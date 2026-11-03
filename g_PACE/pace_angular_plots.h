#ifndef PACE_ANGULAR_PLOTS_H
#define PACE_ANGULAR_PLOTS_H

#include <QString>

#include <array>
#include <vector>

class QWidget;

enum PaceAngularPlotKind
{
    PacePlotEnergyVsTheta = 0,
    PacePlotCountsVsTheta = 1,
    PacePlotCountsVsEnergy = 2,
    PacePlotCrossSectionVsTheta = 3
};

enum PaceSpectrumFitKind
{
    PaceFitNone = 0,
    PaceFitGaussian,
    PaceFitEvaporationNeutron,
    PaceFitEvaporationProton,
    PaceFitEvaporationAlpha
};

struct PaceAngularTablePlotData
{
    QString label;
    double energyLow = 0.0;
    double energyStep = 1.0;
    double angleStep = 1.0;
    double sigmaTotalMb = 0.0;
    int cascades = 1;
    int particleMode = 0;
    bool energyOnly = false;
    double paceCoulombBarrierMeV = 0.0;
    PaceSpectrumFitKind fitKind = PaceFitGaussian;
    std::array<std::array<int, 38>, 33> counts{};
};

struct PaceCMSpectraPlotData
{
    std::array<std::array<int, 52>, 5> counts{};
    std::array<double, 5> averageEnergy{};
    bool valid = false;
};

struct PacePlotData
{
    std::vector<PaceAngularTablePlotData> angularTables;
    PaceCMSpectraPlotData cmSpectra;
};

void paceResetPlotData();

int paceRegisterAngularTable(const QString &label,
                             int tableIndex,
                             const int counts[17][33][38],
                             double energyLow,
                             double energyStep,
                             double angleStep,
                             double sigmaTotalMb,
                             int cascades);

int paceRegisterParticleAngularTable(const QString &label,
                                     int particleMode,
                                     const int counts[4][33][19],
                                     double paceCoulombBarrierMeV,
                                     double sigmaTotalMb,
                                     int cascades);

int paceRegisterGammaEnergyTable(const QString &label,
                                 const int totalGammaCounts[33]);

void paceRegisterCMSpectra(const int counts[5][52],
                           const double averageEnergy[5]);

PacePlotData pacePlotDataSnapshot();

void openPaceAngularPlotWindow(QWidget *parent,
                               const PacePlotData &data,
                               bool plotAllTables,
                               int tableIndex,
                               PaceAngularPlotKind plotKind);

void openPaceCMSpectraPlotWindow(QWidget *parent,
                                 const PacePlotData &data);

#endif // PACE_ANGULAR_PLOTS_H
