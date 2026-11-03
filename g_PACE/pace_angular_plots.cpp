#include "pace_angular_plots.h"

#include <QBoxLayout>
#include <QFileDialog>
#include <QIcon>
#include <QLabel>
#include <QPainter>
#include <QPainterPath>
#include <QPixmap>
#include <QPushButton>
#include <QScrollArea>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

PacePlotData currentPlotData;

constexpr double pi = 3.14159265358979323846;

double tableTotal(const PaceAngularTablePlotData &table)
{
    double total = 0.0;
    for (int angle = 1; angle <= 37; ++angle)
        total += table.counts[32][angle];
    return total;
}

int lastPopulatedAngle(const PaceAngularTablePlotData &table)
{
    for (int angle = 37; angle >= 1; --angle)
        if (table.counts[32][angle] > 0) return angle;
    return 0;
}

double displayedMaxTheta(const PaceAngularTablePlotData &table)
{
    const double rawMax = lastPopulatedAngle(table) * table.angleStep;
    if (rawMax > 36.0) return std::ceil(rawMax);
    if (rawMax > 18.0) return 36.0;
    return 18.0;
}

QString plotKindLabel(PaceAngularPlotKind kind)
{
    if (kind == PacePlotCountsVsTheta) return "N vs θ";
    if (kind == PacePlotCountsVsEnergy) return "N vs Energy";
    if (kind == PacePlotCrossSectionVsTheta) return "dσ/dθ vs θ";
    return "E vs θ";
}

struct GaussianStats
{
    double amplitude = 0.0;
    double mean = 0.0;
    double sigma = 0.0;
    double reducedPoissonDeviance = 0.0;
    double fittedCounts = 0.0;
    QString unavailableReason;
    bool valid = false;
};

struct HistogramStats
{
    double mean = 0.0;
    double median = 0.0;
    double area = 0.0;
    bool valid = false;
};

HistogramStats calculateHistogramStats(const std::vector<double> &x,
                                       const std::vector<double> &values,
                                       double binWidth)
{
    HistogramStats stats;
    if (x.size() != values.size() || x.empty()) return stats;
    double total = 0.0;
    double weighted = 0.0;
    for (int i = 0; i < int(values.size()); ++i)
    {
        const double value = std::max(0.0, values[i]);
        total += value;
        weighted += value * x[i];
    }
    if (total <= 0.0) return stats;
    stats.mean = weighted / total;
    double cumulative = 0.0;
    for (int i = 0; i < int(values.size()); ++i)
    {
        cumulative += std::max(0.0, values[i]);
        if (cumulative >= total * 0.5)
        {
            stats.median = x[i];
            break;
        }
    }
    stats.area = total * std::max(binWidth, 1.0e-12);
    stats.valid = true;
    return stats;
}

double gaussianValue(double amplitude, double mean, double sigma, double x)
{
    if (sigma <= 1.0e-12) return 0.0;
    const double z = (x - mean) / sigma;
    return amplitude * std::exp(-0.5 * z * z);
}

double poissonDevianceTerm(double observed, double expected)
{
    const double mu = std::max(expected, 1.0e-12);
    const double n = std::max(0.0, observed);
    if (n <= 0.0) return 2.0 * mu;
    return 2.0 * (mu - n + n * std::log(n / mu));
}

double gaussianBinProbability(double lower, double upper, double mean, double sigma)
{
    if (sigma <= 1.0e-12 || upper <= lower) return 0.0;
    const double scale = std::sqrt(2.0) * sigma;
    return 0.5 * (std::erf((upper - mean) / scale) -
                  std::erf((lower - mean) / scale));
}

GaussianStats fitGaussianBins(const std::vector<double> &lowerEdges,
                              const std::vector<double> &upperEdges,
                              const std::vector<double> &y,
                              double displayBinWidth)
{
    GaussianStats stats;
    if (lowerEdges.size() != y.size() || upperEdges.size() != y.size() || y.empty())
        return stats;

    int first = -1;
    int last = -1;
    int populatedBins = 0;
    double total = 0.0;
    double weightedX = 0.0;
    double maxY = 0.0;
    for (int i = 0; i < int(y.size()); ++i)
    {
        if (y[i] <= 0.0) continue;
        populatedBins++;
        if (first < 0) first = i;
        last = i;
        total += y[i];
        weightedX += y[i] * 0.5 * (lowerEdges[i] + upperEdges[i]);
        maxY = std::max(maxY, y[i]);
    }
    if (first < 0 || last < first || total <= 0.0)
    {
        stats.unavailableReason = "Gaussian fit unavailable: no counts";
        return stats;
    }
    if (total < 20.0 || populatedBins < 4)
    {
        stats.unavailableReason = "Gaussian fit unavailable: low statistics";
        return stats;
    }

    const double binWidth = std::max(1.0e-9, displayBinWidth);
    const double mean0 = weightedX / total;
    double variance = 0.0;
    for (int i = first; i <= last; ++i)
    {
        const double center = 0.5 * (lowerEdges[i] + upperEdges[i]);
        variance += y[i] * (center - mean0) * (center - mean0);
    }
    const double sigma0 = std::max(binWidth * 0.35, std::sqrt(std::max(variance / total, 1.0e-12)));

    double bestA = maxY;
    double bestMean = mean0;
    double bestSigma = sigma0;
    double bestDeviance = std::numeric_limits<double>::max();

    auto consider = [&](double mean, double sigma)
    {
        if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 1.0e-12) return;
        double shapeSum = 0.0;
        for (int i = 0; i < int(y.size()); ++i)
            shapeSum += gaussianBinProbability(lowerEdges[i], upperEdges[i], mean, sigma);
        if (shapeSum <= 1.0e-18) return;
        const double fittedCounts = total / shapeSum;
        double deviance = 0.0;
        for (int i = 0; i < int(y.size()); ++i)
        {
            const double expected = fittedCounts *
                                    gaussianBinProbability(lowerEdges[i], upperEdges[i], mean, sigma);
            deviance += poissonDevianceTerm(y[i], expected);
        }
        if (deviance < bestDeviance)
        {
            bestDeviance = deviance;
            bestA = fittedCounts * binWidth / (sigma * std::sqrt(2.0 * pi));
            bestMean = mean;
            bestSigma = sigma;
        }
    };

    const double xMin = lowerEdges[first];
    const double xMax = upperEdges[last];
    const double range = std::max(binWidth, xMax - xMin);
    const double sigmaMin = std::max(binWidth * 0.35, 1.0e-9);
    const double sigmaMax = std::max(sigmaMin * 1.05, range * 1.5);
    for (int mi = 0; mi <= 80; ++mi)
    {
        const double mean = xMin + (xMax - xMin) * double(mi) / 80.0;
        for (int si = 0; si <= 80; ++si)
        {
            const double frac = double(si) / 80.0;
            consider(mean, sigmaMin * std::pow(sigmaMax / sigmaMin, frac));
        }
    }

    for (int pass = 0; pass < 3; ++pass)
    {
        const double meanWindow = std::max(0.25 * binWidth, range / std::pow(4.0, pass + 1.0));
        const double sigmaScale = std::pow(1.8, 1.0 / (pass + 1.0));
        const double localMeanMin = std::max(xMin - binWidth, bestMean - meanWindow);
        const double localMeanMax = std::min(xMax + binWidth, bestMean + meanWindow);
        const double localSigmaMin = std::max(sigmaMin, bestSigma / sigmaScale);
        const double localSigmaMax = std::min(sigmaMax, bestSigma * sigmaScale);
        for (int mi = 0; mi <= 40; ++mi)
        {
            const double mean = localMeanMin + (localMeanMax - localMeanMin) * double(mi) / 40.0;
            for (int si = 0; si <= 40; ++si)
            {
                const double frac = double(si) / 40.0;
                consider(mean, localSigmaMin *
                                   std::pow(std::max(localSigmaMax / localSigmaMin, 1.0), frac));
            }
        }
    }

    stats.amplitude = bestA;
    stats.mean = bestMean;
    stats.sigma = bestSigma;
    stats.fittedCounts = bestA * bestSigma * std::sqrt(2.0 * pi) / binWidth;
    stats.reducedPoissonDeviance = bestDeviance / std::max(1, int(y.size()) - 3);
    stats.valid = std::isfinite(stats.amplitude) && std::isfinite(stats.mean) &&
                  std::isfinite(stats.sigma) && stats.amplitude > 0.0 && stats.sigma > 0.0;
    if (!stats.valid) stats.unavailableReason = "Gaussian fit unavailable: numerical failure";
    return stats;
}

QString compactNumber(double value)
{
    if (!std::isfinite(value)) return "nan";
    if (std::fabs(value) >= 1000.0 || (std::fabs(value) > 0.0 && std::fabs(value) < 0.01))
        return QString::number(value, 'e', 2);
    return QString::number(value, 'g', 4);
}

struct EvaporationStats
{
    bool valid = false;
    QString unavailableReason;
    double normalization = 0.0;
    double temperatureMeV = 0.0;
    double barrierMeV = 0.0;
    double reducedPoissonDeviance = 0.0;
    double reducedPearsonChiSquare = 0.0;
};

bool isChargedFit(PaceSpectrumFitKind kind)
{
    return kind == PaceFitEvaporationProton || kind == PaceFitEvaporationAlpha;
}

double evaporationDensity(double energy,
                          PaceSpectrumFitKind kind,
                          double temperature,
                          double barrier)
{
    if (energy < 0.0 || temperature <= 0.0) return 0.0;
    double penetrability = 1.0;
    if (isChargedFit(kind))
    {
        const double diffuseness = kind == PaceFitEvaporationAlpha ? 1.5 : 1.0;
        const double exponent = std::clamp((barrier - energy) / diffuseness, -700.0, 700.0);
        penetrability = 1.0 / (1.0 + std::exp(exponent));
    }
    return std::sqrt(std::max(energy, 0.0)) * penetrability *
           std::exp(std::clamp(-energy / temperature, -700.0, 0.0));
}

EvaporationStats fitEvaporationSpectrum(const std::vector<double> &x,
                                        const std::vector<double> &y,
                                        PaceSpectrumFitKind kind,
                                        double paceBarrier)
{
    EvaporationStats result;
    if (x.size() != y.size() || x.empty() || kind < PaceFitEvaporationNeutron)
    {
        result.unavailableReason = "Evaporation fit unavailable";
        return result;
    }

    double total = 0.0;
    int populated = 0;
    int first = -1;
    int last = -1;
    for (int i = 0; i < int(y.size()); ++i)
    {
        total += std::max(0.0, y[i]);
        if (y[i] > 0.0)
        {
            populated++;
            if (first < 0) first = i;
            last = i;
        }
    }
    if (total < 20.0 || populated < 4 || first < 0)
    {
        result.unavailableReason = "Evaporation fit unavailable: low statistics";
        return result;
    }

    if (total >= 50.0 && populated >= 6)
    {
        double cumulative = 0.0;
        int coreFirst = first;
        int coreLast = last;
        for (int i = 0; i < int(y.size()); ++i)
        {
            cumulative += std::max(0.0, y[i]);
            if (cumulative >= 0.005 * total) { coreFirst = i; break; }
        }
        cumulative = 0.0;
        for (int i = 0; i < int(y.size()); ++i)
        {
            cumulative += std::max(0.0, y[i]);
            if (cumulative >= 0.995 * total) { coreLast = i; break; }
        }
        first = std::max(0, coreFirst - 1);
        last = std::min(int(y.size()) - 1, coreLast + 1);
    }

    const bool charged = isChargedFit(kind);
    if (charged && paceBarrier <= 0.0)
        paceBarrier = kind == PaceFitEvaporationAlpha ? 10.0 : 5.0;
    const double barrierMin = charged ? std::max(0.0, std::min(0.5 * paceBarrier, paceBarrier - 4.0)) : 0.0;
    const double barrierMax = charged
                                  ? std::max(barrierMin + 0.5,
                                             std::min(kind == PaceFitEvaporationAlpha ? 60.0 : 30.0,
                                                      std::max(1.6 * paceBarrier, paceBarrier + 4.0)))
                                  : 0.0;

    double bestT = 0.0;
    double bestBarrier = 0.0;
    double bestNormalization = 0.0;
    double bestDeviance = std::numeric_limits<double>::max();
    auto consider = [&](double temperature, double barrier)
    {
        double modelSum = 0.0;
        double observedSum = 0.0;
        for (int i = first; i <= last; ++i)
        {
            modelSum += evaporationDensity(x[i], kind, temperature, barrier);
            observedSum += std::max(0.0, y[i]);
        }
        if (modelSum <= 1.0e-18 || observedSum <= 0.0) return;
        const double normalization = observedSum / modelSum;
        double deviance = 0.0;
        for (int i = first; i <= last; ++i)
            deviance += poissonDevianceTerm(y[i],
                                            normalization * evaporationDensity(x[i], kind, temperature, barrier));
        if (deviance < bestDeviance)
        {
            bestDeviance = deviance;
            bestT = temperature;
            bestBarrier = barrier;
            bestNormalization = normalization;
        }
    };

    for (int ti = 0; ti <= 240; ++ti)
    {
        const double temperature = 0.2 * std::pow(10.0 / 0.2, double(ti) / 240.0);
        const int barrierSteps = charged ? 80 : 0;
        for (int bi = 0; bi <= barrierSteps; ++bi)
        {
            const double barrier = charged
                                       ? barrierMin + (barrierMax - barrierMin) * double(bi) / barrierSteps
                                       : 0.0;
            consider(temperature, barrier);
        }
    }

    for (int pass = 0; pass < 4 && bestNormalization > 0.0; ++pass)
    {
        const double tWindow = std::max(0.05, bestT * 0.18 * std::pow(0.35, pass));
        const double bWindow = charged
                                   ? std::max(0.05, (barrierMax - barrierMin) * 0.12 * std::pow(0.35, pass))
                                   : 0.0;
        const double tLow = std::max(0.2, bestT - tWindow);
        const double tHigh = std::min(10.0, bestT + tWindow);
        const double bLow = charged ? std::max(barrierMin, bestBarrier - bWindow) : 0.0;
        const double bHigh = charged ? std::min(barrierMax, bestBarrier + bWindow) : 0.0;
        for (int ti = 0; ti <= 80; ++ti)
            for (int bi = 0; bi <= (charged ? 40 : 0); ++bi)
                consider(tLow + (tHigh - tLow) * double(ti) / 80.0,
                         charged ? bLow + (bHigh - bLow) * double(bi) / 40.0 : 0.0);
    }

    if (bestNormalization <= 0.0 || !std::isfinite(bestDeviance))
    {
        result.unavailableReason = "Evaporation fit unavailable";
        return result;
    }
    const int freeParameters = charged ? 3 : 2;
    result.valid = true;
    result.normalization = bestNormalization;
    result.temperatureMeV = bestT;
    result.barrierMeV = bestBarrier;
    result.reducedPoissonDeviance = bestDeviance / std::max(1, last - first + 1 - freeParameters);
    double pearsonChiSquare = 0.0;
    for (int i = first; i <= last; ++i)
    {
        const double expected = std::max(1.0e-12,
                                         bestNormalization * evaporationDensity(
                                             x[i], kind, bestT, bestBarrier));
        const double difference = y[i] - expected;
        pearsonChiSquare += difference * difference / std::max(1.0, expected);
    }
    result.reducedPearsonChiSquare =
        pearsonChiSquare / std::max(1, last - first + 1 - freeParameters);
    return result;
}

class PaceEnergyThetaPlot : public QWidget
{
public:
    explicit PaceEnergyThetaPlot(const PaceAngularTablePlotData &table, QWidget *parent = nullptr)
        : QWidget(parent), table(table)
    {
        setMinimumSize(980, 560);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing, true);
        painter.fillRect(rect(), Qt::white);

        const QRect plotRect(90, 35, width() - 210, height() - 115);
        painter.setPen(Qt::black);
        painter.drawRect(plotRect);
        if (tableTotal(table) <= 0.0)
        {
            painter.drawText(plotRect, Qt::AlignCenter, "No histogram data");
            return;
        }

        const double maxTheta = displayedMaxTheta(table);
        const int thetaBins = std::max(1, std::min(37, int(std::ceil(maxTheta / table.angleStep))));
        int highestEnergyBin = 1;
        int maxCount = 0;
        for (int energy = 1; energy <= 31; ++energy)
            for (int angle = 1; angle <= thetaBins; ++angle)
            {
                if (table.counts[energy][angle] > 0) highestEnergyBin = std::max(highestEnergyBin, energy);
                maxCount = std::max(maxCount, table.counts[energy][angle]);
            }

        const double cellWidth = double(plotRect.width()) / thetaBins;
        const double cellHeight = double(plotRect.height()) / highestEnergyBin;
        for (int energy = 1; energy <= highestEnergyBin; ++energy)
            for (int angle = 1; angle <= thetaBins; ++angle)
            {
                const int count = table.counts[energy][angle];
                const double fraction = maxCount > 0 ? double(count) / maxCount : 0.0;
                QColor color = count > 0
                                   ? QColor(int(255 * fraction), int(230 * (1.0 - fraction)), 255 - int(220 * fraction))
                                   : QColor(250, 250, 250);
                const int x = plotRect.left() + int((angle - 1) * cellWidth);
                const int y = plotRect.bottom() - int(energy * cellHeight);
                QRect cell(x, y, int(cellWidth) + 1, int(cellHeight) + 1);
                painter.fillRect(cell, color);
                painter.setPen(QColor(220, 220, 220));
                painter.drawRect(cell);
                if (count > 0)
                {
                    painter.setPen(Qt::black);
                    painter.drawText(cell, Qt::AlignCenter, QString::number(count));
                }
            }

        painter.setPen(QPen(Qt::black, 2));
        painter.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        painter.drawLine(plotRect.bottomLeft(), plotRect.topLeft());
        const int thetaTick = maxTheta > 36.0 ? 20 : (maxTheta > 18.0 ? 5 : 2);
        for (int theta = 0; theta <= int(maxTheta); theta += thetaTick)
        {
            const int x = plotRect.left() + int(theta / maxTheta * plotRect.width());
            painter.drawLine(x, plotRect.bottom(), x, plotRect.bottom() + 6);
            painter.drawText(x - 20, plotRect.bottom() + 8, 40, 18, Qt::AlignCenter, QString::number(theta));
        }
        for (int energy = 0; energy <= highestEnergyBin; ++energy)
        {
            if (energy != 0 && energy != 1 && energy != highestEnergyBin && (energy - 1) % 5 != 0) continue;
            const int y = plotRect.bottom() - int(double(energy) / highestEnergyBin * plotRect.height());
            double value = energy == 0 ? 0.0 : (energy == 1 ? table.energyLow : table.energyLow + (energy - 1) * table.energyStep);
            painter.drawLine(plotRect.left() - 6, y, plotRect.left(), y);
            painter.drawText(5, y - 10, 75, 20, Qt::AlignRight | Qt::AlignVCenter, QString::number(value, 'f', 1));
        }

        painter.drawText(plotRect.left(), height() - 35, plotRect.width(), 24, Qt::AlignCenter, "θ range (deg)");
        painter.save();
        painter.translate(28, plotRect.center().y());
        painter.rotate(-90);
        painter.drawText(QRect(-plotRect.height() / 2, -20, plotRect.height(), 20), Qt::AlignCenter, "Energy range (MeV)");
        painter.restore();

        const int legendX = plotRect.right() + 25;
        const int legendY = plotRect.top() + 20;
        for (int i = 0; i < 180; ++i)
        {
            const double fraction = 1.0 - double(i) / 179.0;
            painter.setPen(QColor(int(255 * fraction), int(230 * (1.0 - fraction)), 255 - int(220 * fraction)));
            painter.drawLine(legendX, legendY + i, legendX + 22, legendY + i);
        }
        painter.setPen(Qt::black);
        painter.drawRect(legendX, legendY, 22, 180);
        painter.drawText(legendX + 30, legendY - 4, 70, 20, Qt::AlignLeft, QString::number(maxCount));
        painter.drawText(legendX + 30, legendY + 165, 70, 20, Qt::AlignLeft, "0");
        painter.drawText(legendX - 8, legendY + 195, 90, 20, Qt::AlignCenter, "Counts");
    }

private:
    PaceAngularTablePlotData table;
};

class PaceOneDimensionalPlot : public QWidget
{
public:
    PaceOneDimensionalPlot(const PaceAngularTablePlotData &table,
                           PaceAngularPlotKind kind,
                           QWidget *parent = nullptr)
        : QWidget(parent), table(table), kind(kind)
    {
        setMinimumSize(980, 560);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing, true);
        painter.fillRect(rect(), Qt::white);
        const QRect plotRect(90, 45, width() - 150, height() - 130);
        painter.setPen(Qt::black);
        painter.drawRect(plotRect);

        std::vector<double> x;
        std::vector<double> values;
        std::vector<double> binLower;
        std::vector<double> binUpper;
        double minX = 0.0;
        double maxX = 1.0;
        double histogramBinWidth = 1.0;
        double separatedOverflowCount = 0.0;
        bool alphaOverflowSeparated = false;
        bool displayedEnergyIncludesOverflow = false;
        QString title;
        QString xTitle;
        QString yTitle;

        if (kind == PacePlotCountsVsEnergy)
        {
            if (table.particleMode == 0)
            {
                // Residual bin 1 is the aggregate "Below ELOW" count and bin
                // 31 is the aggregate "Above" count. They have no finite bin
                // width, so they must not be drawn or included in a Gaussian fit.
                constexpr int firstFiniteBin = 2;
                constexpr int lastFiniteBin = 30;
                constexpr int finiteBinCount = lastFiniteBin - firstFiniteBin + 1;
                x.resize(finiteBinCount);
                values.assign(finiteBinCount, 0.0);
                binLower.resize(finiteBinCount);
                binUpper.resize(finiteBinCount);
                for (int energy = firstFiniteBin; energy <= lastFiniteBin; ++energy)
                {
                    const int outputBin = energy - firstFiniteBin;
                    for (int angle = 1; angle <= 37; ++angle)
                        values[outputBin] += table.counts[energy][angle];
                    binLower[outputBin] = table.energyLow + (energy - 2.0) * table.energyStep;
                    binUpper[outputBin] = binLower[outputBin] + table.energyStep;
                    x[outputBin] = 0.5 * (binLower[outputBin] + binUpper[outputBin]);
                }
                minX = binLower.front();
                maxX = std::max(binUpper.back(), minX + table.energyStep);
            }
            else
            {
                x.resize(31);
                values.assign(31, 0.0);
                binLower.resize(31);
                binUpper.resize(31);
                displayedEnergyIncludesOverflow = true;
                for (int energy = 1; energy <= 31; ++energy)
                {
                    for (int angle = 1; angle <= 37; ++angle)
                        values[energy - 1] += table.counts[energy][angle];
                    if (energy == 1)
                    {
                        binLower[energy - 1] = 0.0;
                        binUpper[energy - 1] = table.energyLow;
                    }
                    else
                    {
                        binLower[energy - 1] = table.energyLow + (energy - 2.0) * table.energyStep;
                        binUpper[energy - 1] = binLower[energy - 1] + table.energyStep;
                    }
                    x[energy - 1] = 0.5 * (binLower[energy - 1] + binUpper[energy - 1]);
                }
                maxX = std::max(binUpper.back(), 1.0);
            }
            histogramBinWidth = table.energyStep;
            xTitle = table.particleMode > 0 ? "Laboratory energy (MeV)" : "Energy (MeV)";
            if (table.fitKind == PaceFitEvaporationAlpha && !values.empty())
            {
                alphaOverflowSeparated = true;
                separatedOverflowCount = values.back();
                values.pop_back();
                x.pop_back();
                binLower.pop_back();
                binUpper.pop_back();
                displayedEnergyIncludesOverflow = false;
                maxX = std::max(binUpper.back(), 1.0);
                xTitle = "Lab energy (MeV)";
            }
            title = "N = f(Energy)";
            yTitle = "Counts N";
        }
        else
        {
            maxX = displayedMaxTheta(table);
            const int bins = std::max(1, std::min(37, int(std::ceil(maxX / table.angleStep))));
            histogramBinWidth = table.angleStep;
            x.resize(bins);
            values.assign(bins, 0.0);
            binLower.resize(bins);
            binUpper.resize(bins);
            for (int angle = 1; angle <= bins; ++angle)
            {
                binLower[angle - 1] = (angle - 1.0) * table.angleStep;
                binUpper[angle - 1] = angle * table.angleStep;
                x[angle - 1] = 0.5 * (binLower[angle - 1] + binUpper[angle - 1]);
                const double count = table.counts[32][angle];
                values[angle - 1] = kind == PacePlotCrossSectionVsTheta
                                        ? table.sigmaTotalMb * count /
                                              (std::max(1, table.cascades) * table.angleStep)
                                        : count;
            }
            title = kind == PacePlotCrossSectionVsTheta ? "dσ/dθ = f(θ)" : "N = f(θ)";
            xTitle = "θ (deg)";
            yTitle = kind == PacePlotCrossSectionVsTheta ? "dσ/dθ (mb/deg)" : "Counts N";
        }

        const double maxYValue = values.empty() ? 0.0 : *std::max_element(values.begin(), values.end());
        GaussianStats gaussian;
        EvaporationStats evaporation;
        std::vector<double> fitX = x;
        std::vector<double> fitValues = values;
        std::vector<double> fitLower = binLower;
        std::vector<double> fitUpper = binUpper;
        // A displayed particle overflow bin is open-ended. Show it, but never
        // assign its count to an arbitrary finite center during a fit.
        if (kind == PacePlotCountsVsEnergy && displayedEnergyIncludesOverflow && fitX.size() > 1)
        {
            fitX.pop_back();
            fitValues.pop_back();
            fitLower.pop_back();
            fitUpper.pop_back();
        }
        if (kind == PacePlotCountsVsEnergy && table.fitKind == PaceFitGaussian)
            gaussian = fitGaussianBins(fitLower, fitUpper, fitValues, histogramBinWidth);
        else if (kind == PacePlotCountsVsEnergy && table.fitKind >= PaceFitEvaporationNeutron)
            evaporation = fitEvaporationSpectrum(fitX, fitValues, table.fitKind,
                                                 table.paceCoulombBarrierMeV);
        const HistogramStats histogram = calculateHistogramStats(x, values, histogramBinWidth);
        double maxY = std::max(1.0, maxYValue * 1.15);
        if (gaussian.valid) maxY = std::max(maxY, gaussian.amplitude * 1.15);
        if (evaporation.valid)
            for (double energy = minX; energy <= maxX; energy += 0.1)
                maxY = std::max(maxY, evaporation.normalization *
                                          evaporationDensity(energy, table.fitKind,
                                                             evaporation.temperatureMeV,
                                                             evaporation.barrierMeV) * 1.15);

        auto mapX = [&](double value) {
            return plotRect.left() + int((value - minX) / std::max(1.0e-12, maxX - minX) * plotRect.width());
        };
        auto mapY = [&](double value) {
            return plotRect.bottom() - int(value / maxY * plotRect.height());
        };

        painter.setPen(QPen(QColor(232, 232, 232), 1));
        for (int i = 0; i <= 5; ++i)
        {
            const int y = plotRect.bottom() - i * plotRect.height() / 5;
            painter.drawLine(plotRect.left(), y, plotRect.right(), y);
            painter.setPen(Qt::black);
            painter.drawText(5, y - 9, 75, 18, Qt::AlignRight | Qt::AlignVCenter,
                             compactNumber(maxY * i / 5.0));
            painter.setPen(QPen(QColor(232, 232, 232), 1));
        }

        for (int i = 0; i < int(values.size()); ++i)
        {
            const int left = mapX(binLower[i]);
            const int right = mapX(binUpper[i]);
            const int top = mapY(values[i]);
            QRect bar(left + 1, top, std::max(1, right - left - 2), std::max(0, plotRect.bottom() - top));
            painter.fillRect(bar, values[i] > 0.0 ? QColor(160, 190, 255) : QColor(248, 248, 248));
            painter.setPen(Qt::black);
            painter.drawRect(bar);
        }

        if (gaussian.valid)
        {
            QPainterPath curve;
            for (int i = 0; i <= 500; ++i)
            {
                const double value = minX + (maxX - minX) * i / 500.0;
                const QPointF point(mapX(value), mapY(gaussianValue(gaussian.amplitude, gaussian.mean, gaussian.sigma, value)));
                if (i == 0) curve.moveTo(point); else curve.lineTo(point);
            }
            painter.setPen(QPen(QColor(220, 40, 40, 155), 2, Qt::DotLine));
            painter.drawPath(curve);
        }
        if (evaporation.valid)
        {
            QPainterPath curve;
            for (int i = 0; i <= 500; ++i)
            {
                const double energy = minX + (maxX - minX) * i / 500.0;
                const double expected = evaporation.normalization *
                                        evaporationDensity(energy, table.fitKind,
                                                           evaporation.temperatureMeV,
                                                           evaporation.barrierMeV);
                const QPointF point(mapX(energy), mapY(expected));
                if (i == 0) curve.moveTo(point); else curve.lineTo(point);
            }
            painter.setPen(QPen(QColor(210, 75, 30), 2));
            painter.drawPath(curve);
            if (isChargedFit(table.fitKind) &&
                (!alphaOverflowSeparated ||
                 (evaporation.barrierMeV >= minX && evaporation.barrierMeV <= maxX)))
            {
                painter.setPen(QPen(QColor(80, 80, 80), 2, Qt::DashLine));
                const int barrierX = mapX(evaporation.barrierMeV);
                painter.drawLine(barrierX, plotRect.top(), barrierX, plotRect.bottom());
            }
        }

        painter.setPen(QPen(Qt::black, 2));
        painter.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());
        painter.drawLine(plotRect.bottomLeft(), plotRect.topLeft());
        for (int i = 0; i <= 10; ++i)
        {
            const double value = minX + (maxX - minX) * i / 10.0;
            const int px = mapX(value);
            painter.drawLine(px, plotRect.bottom(), px, plotRect.bottom() + 5);
            painter.drawText(px - 25, plotRect.bottom() + 7, 50, 18, Qt::AlignCenter,
                             QString::number(value, 'g', 4));
        }

        QFont titleFont = painter.font();
        titleFont.setBold(true);
        titleFont.setPointSize(11);
        painter.setFont(titleFont);
        painter.drawText(plotRect.left(), 8, plotRect.width(), 24, Qt::AlignCenter, title);
        titleFont.setBold(false);
        titleFont.setPointSize(9);
        painter.setFont(titleFont);
        painter.drawText(plotRect.left(), height() - 34, plotRect.width(), 22, Qt::AlignCenter, xTitle);
        painter.save();
        painter.translate(28, plotRect.center().y());
        painter.rotate(-90);
        painter.drawText(QRect(-plotRect.height() / 2, -20, plotRect.height(), 20), Qt::AlignCenter, yTitle);
        painter.restore();

        const bool fitUnavailable = kind == PacePlotCountsVsEnergy &&
                                    ((!gaussian.unavailableReason.isEmpty() && table.fitKind == PaceFitGaussian) ||
                                     (!evaporation.unavailableReason.isEmpty() &&
                                      table.fitKind >= PaceFitEvaporationNeutron));
        const bool showHistogramStats = histogram.valid &&
                                        (kind == PacePlotCrossSectionVsTheta ||
                                         (kind == PacePlotCountsVsEnergy &&
                                          !gaussian.valid && !evaporation.valid));
        const int fitHeight = gaussian.valid ? 105
                                             : (evaporation.valid
                                                    ? (alphaOverflowSeparated ? 105 : 90)
                                                    : (fitUnavailable ? 35 : 0));
        const int legendHeight = 36 + (showHistogramStats ? 45 : 0) + fitHeight;
        const int legendWidth = evaporation.valid ? 300 : 230;
        const QRect legend(plotRect.right() - legendWidth - 12, plotRect.top() + 12,
                           legendWidth, legendHeight);
        painter.fillRect(legend, QColor(255, 255, 255, 225));
        painter.setPen(QColor(120, 120, 120));
        painter.drawRect(legend);
        QFont legendFont = painter.font();
        legendFont.setPointSize(8);
        legendFont.setBold(false);
        painter.setFont(legendFont);
        painter.setPen(Qt::black);
        const double legendCount = kind == PacePlotCountsVsEnergy && table.particleMode == 0
                                       ? tableTotal(table)
                                       : (histogram.valid ? histogram.area / histogramBinWidth
                                                          : tableTotal(table)) +
                                             separatedOverflowCount;
        painter.drawText(legend.left() + 10, legend.top() + 8, legend.width() - 20, 14,
                         Qt::AlignLeft, "N = " + compactNumber(legendCount));
        int cursorY = legend.top() + 33;
        if (showHistogramStats)
        {
            painter.drawText(legend.left() + 10, cursorY, legend.width() - 20, 14, Qt::AlignLeft,
                             "Mean = " + compactNumber(histogram.mean));
            painter.drawText(legend.left() + 10, cursorY + 15, legend.width() - 20, 14, Qt::AlignLeft,
                             "Median = " + compactNumber(histogram.median));
            painter.drawText(legend.left() + 10, cursorY + 30, legend.width() - 20, 14, Qt::AlignLeft,
                             "Area = " + compactNumber(histogram.area));
            cursorY += 45;
        }
        if (gaussian.valid)
        {
            painter.setPen(QPen(QColor(220, 40, 40, 155), 2, Qt::DotLine));
            painter.drawLine(legend.left() + 10, cursorY, legend.left() + 42, cursorY);
            painter.setPen(Qt::black);
            painter.drawText(legend.left() + 50, cursorY - 9, 150, 16, Qt::AlignLeft, "Gaussian fit");
            painter.drawText(legend.left() + 10, cursorY + 15, 190, 14, Qt::AlignLeft,
                             "Amp = " + compactNumber(gaussian.amplitude));
            painter.drawText(legend.left() + 10, cursorY + 30, 190, 14, Qt::AlignLeft,
                             "Mean = " + compactNumber(gaussian.mean));
            painter.drawText(legend.left() + 10, cursorY + 45, 190, 14, Qt::AlignLeft,
                             "Median = " + compactNumber(histogram.median));
            painter.drawText(legend.left() + 10, cursorY + 60, 190, 14, Qt::AlignLeft,
                             "Poisson D/DoF = " + compactNumber(gaussian.reducedPoissonDeviance));
            painter.drawText(legend.left() + 10, cursorY + 75, 190, 14, Qt::AlignLeft,
                             "Fit counts = " + compactNumber(gaussian.fittedCounts));
        }
        else if (evaporation.valid)
        {
            painter.setPen(QPen(QColor(210, 75, 30), 2));
            painter.drawLine(legend.left() + 10, cursorY, legend.left() + 42, cursorY);
            painter.setPen(Qt::black);
            painter.drawText(legend.left() + 50, cursorY - 9, legend.width() - 60, 16,
                             Qt::AlignLeft,
                             alphaOverflowSeparated ? "Boltzmann fit" : "Evaporation fit");
            int infoY = cursorY + 15;
            if (alphaOverflowSeparated)
            {
                painter.drawText(legend.left() + 10, infoY, legend.width() - 20, 14,
                                 Qt::AlignLeft,
                                 "A = " + compactNumber(evaporation.normalization));
                infoY += 15;
            }
            painter.drawText(legend.left() + 10, infoY, legend.width() - 20, 14,
                             Qt::AlignLeft, "T = " + compactNumber(evaporation.temperatureMeV) + " MeV");
            infoY += 15;
            if (isChargedFit(table.fitKind))
            {
                painter.drawText(legend.left() + 10, infoY, legend.width() - 20, 14,
                                 Qt::AlignLeft, "B = " +
                                                   compactNumber(evaporation.barrierMeV) + " MeV");
                infoY += 15;
            }
            painter.drawText(legend.left() + 10, infoY, legend.width() - 20, 14,
                             Qt::AlignLeft,
                             alphaOverflowSeparated
                                 ? "Chi^2(red) = " +
                                       compactNumber(evaporation.reducedPearsonChiSquare)
                                 : "Poisson D/DoF = " +
                                       compactNumber(evaporation.reducedPoissonDeviance));
        }
        else if (fitUnavailable)
        {
            const QString reason = table.fitKind == PaceFitGaussian
                                       ? gaussian.unavailableReason
                                       : evaporation.unavailableReason;
            painter.drawText(legend.left() + 10, cursorY, legend.width() - 20, 28,
                             Qt::AlignLeft | Qt::TextWordWrap, reason);
        }
    }

private:
    PaceAngularTablePlotData table;
    PaceAngularPlotKind kind;
};

struct CMSpectrumSeries
{
    QString label;
    QColor color;
    std::array<int, 52> counts{};
    double averageEnergy = 0.0;
};

class PaceCMSpectraPlot : public QWidget
{
public:
    explicit PaceCMSpectraPlot(const PaceCMSpectraPlotData &data, QWidget *parent = nullptr)
        : QWidget(parent)
    {
        const std::array<QString, 4> labels = {"Neutrons", "Protons", "Alpha particles", "Gamma particles"};
        const std::array<QColor, 4> colors = {QColor(54, 112, 181), QColor(205, 73, 73),
                                              QColor(65, 150, 90), QColor(176, 125, 38)};
        for (int i = 0; i < 4; ++i)
        {
            CMSpectrumSeries series;
            series.label = labels[i];
            series.color = colors[i];
            series.counts = data.counts[i + 1];
            series.averageEnergy = data.averageEnergy[i + 1];
            spectra.push_back(series);
        }
        setMinimumSize(1120, 760);
    }

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing, true);
        painter.fillRect(rect(), Qt::white);
        QFont font = painter.font();
        font.setPointSize(12);
        font.setBold(true);
        painter.setFont(font);
        painter.drawText(0, 10, width(), 26, Qt::AlignCenter, "C.M. Spectra of Emitted Particles");

        const int margin = 18;
        const int gap = 20;
        const int top = 48;
        const int panelWidth = (width() - 2 * margin - gap) / 2;
        const int panelHeight = (height() - top - margin - gap) / 2;
        for (int i = 0; i < int(spectra.size()); ++i)
        {
            const QRect panel(margin + (i % 2) * (panelWidth + gap),
                              top + (i / 2) * (panelHeight + gap), panelWidth, panelHeight);
            drawPanel(painter, panel, spectra[i]);
        }
    }

private:
    void drawPanel(QPainter &painter, const QRect &panel, const CMSpectrumSeries &series)
    {
        const QRect plot(panel.left() + 58, panel.top() + 34, panel.width() - 72, panel.height() - 92);
        int maxCount = 0;
        for (int bin = 1; bin <= 50; ++bin) maxCount = std::max(maxCount, series.counts[bin]);
        const double maxY = std::max(1.0, std::ceil(maxCount * 1.15));
        auto mapY = [&](double value) { return plot.bottom() - int(value / maxY * plot.height()); };

        QFont font = painter.font();
        font.setPointSize(10);
        font.setBold(true);
        painter.setFont(font);
        painter.setPen(Qt::black);
        painter.drawText(panel.left(), panel.top() + 4, panel.width(), 18, Qt::AlignCenter, series.label);
        font.setPointSize(8);
        font.setBold(false);
        painter.setFont(font);

        for (int i = 0; i <= 4; ++i)
        {
            const int y = plot.bottom() - i * plot.height() / 4;
            painter.setPen(QColor(232, 232, 232));
            painter.drawLine(plot.left(), y, plot.right(), y);
            painter.setPen(Qt::black);
            painter.drawText(panel.left() + 2, y - 9, 48, 18, Qt::AlignRight | Qt::AlignVCenter,
                             QString::number(int(std::round(maxY * i / 4.0))));
        }

        const double barWidth = double(plot.width()) / 50.0;
        QColor fill = series.color;
        fill.setAlpha(170);
        const QColor edge = series.color.darker(115);
        for (int bin = 1; bin <= 50; ++bin)
        {
            const int left = plot.left() + int((bin - 1) * barWidth);
            const int right = plot.left() + int(bin * barWidth);
            const int y = mapY(series.counts[bin]);
            QRect bar(left + 1, y, std::max(1, right - left - 2), std::max(0, plot.bottom() - y));
            if (series.counts[bin] > 0)
            {
                painter.fillRect(bar, fill);
                painter.setPen(edge);
                painter.drawRect(bar);
            }
        }

        painter.setPen(Qt::black);
        painter.drawRect(plot);
        for (int energy = 0; energy <= 50; energy += 10)
        {
            const int x = plot.left() + int(double(energy) / 50.0 * plot.width());
            painter.drawLine(x, plot.bottom(), x, plot.bottom() + 4);
            painter.drawText(x - 18, plot.bottom() + 6, 36, 16, Qt::AlignCenter, QString::number(energy));
        }

        const QRect legend(plot.right() - 148, plot.top() + 8, 140, 50);
        painter.fillRect(legend, QColor(255, 255, 255, 225));
        painter.setPen(QColor(150, 150, 150));
        painter.drawRect(legend);
        painter.fillRect(QRect(legend.left() + 8, legend.top() + 10, 18, 10), fill);
        painter.setPen(Qt::black);
        painter.drawText(legend.left() + 32, legend.top() + 5, 100, 16, Qt::AlignLeft, "Counts/bin");
        painter.drawText(legend.left() + 8, legend.top() + 25, 125, 16, Qt::AlignLeft,
                         "Total " + QString::number(series.counts[51]) +
                             ", Avg " + QString::number(series.averageEnergy, 'f', 2));
        if (series.counts[51] == 0) painter.drawText(plot, Qt::AlignCenter, "No counts");

        painter.drawText(plot.left(), panel.bottom() - 20, plot.width(), 16, Qt::AlignCenter, "Energy bin Ex (MeV)");
        painter.save();
        painter.translate(panel.left() + 16, plot.center().y());
        painter.rotate(-90);
        painter.drawText(QRect(-plot.height() / 2, -20, plot.height(), 20), Qt::AlignCenter, "Counts");
        painter.restore();
    }

    std::vector<CMSpectrumSeries> spectra;
};

QPushButton *addSaveButton(QHBoxLayout *row, QWidget *owner, QWidget *renderTarget)
{
    QPushButton *button = new QPushButton("Save as PNG");
    button->setFixedHeight(26);
    button->setMaximumWidth(110);
    row->addWidget(button, 0, Qt::AlignTop | Qt::AlignRight);
    QObject::connect(button, &QPushButton::clicked, owner, [owner, renderTarget]() {
        const QString fileName = QFileDialog::getSaveFileName(owner, QObject::tr("Save Plot as PNG"),
                                                              QString(), QObject::tr("PNG Image (*.png)"));
        if (fileName.isEmpty()) return;
        QPixmap pixmap(renderTarget->size());
        pixmap.fill(Qt::white);
        renderTarget->render(&pixmap);
        pixmap.save(fileName, "PNG");
    });
    return button;
}

} // namespace

void paceResetPlotData()
{
    currentPlotData = PacePlotData();
}

int paceRegisterAngularTable(const QString &label,
                             int tableIndex,
                             const int counts[17][33][38],
                             double energyLow,
                             double energyStep,
                             double angleStep,
                             double sigmaTotalMb,
                             int cascades)
{
    if (tableIndex < 0 || tableIndex >= 17) return -1;
    PaceAngularTablePlotData table;
    table.label = label;
    table.energyLow = energyLow;
    table.energyStep = energyStep;
    table.angleStep = angleStep;
    table.sigmaTotalMb = sigmaTotalMb;
    table.cascades = std::max(1, cascades);
    for (int energy = 0; energy < 33; ++energy)
        for (int angle = 0; angle < 38; ++angle)
            table.counts[energy][angle] = counts[tableIndex][energy][angle];
    currentPlotData.angularTables.push_back(table);
    return int(currentPlotData.angularTables.size()) - 1;
}

int paceRegisterParticleAngularTable(const QString &label,
                                     int particleMode,
                                     const int counts[4][33][19],
                                     double paceCoulombBarrierMeV,
                                     double sigmaTotalMb,
                                     int cascades)
{
    if (particleMode < 1 || particleMode > 3) return -1;
    PaceAngularTablePlotData table;
    table.label = label;
    // PACE particle bins are 0-1, 1-2, ..., 29-30 MeV plus an overflow bin.
    table.energyLow = 1.0;
    table.energyStep = 1.0;
    table.angleStep = 10.0;
    table.sigmaTotalMb = sigmaTotalMb;
    table.cascades = std::max(1, cascades);
    table.particleMode = particleMode;
    table.paceCoulombBarrierMeV = paceCoulombBarrierMeV;
    table.fitKind = particleMode == 1 ? PaceFitEvaporationNeutron
                                     : (particleMode == 2 ? PaceFitEvaporationProton
                                                          : PaceFitEvaporationAlpha);
    for (int energy = 0; energy < 33; ++energy)
        for (int angle = 0; angle < 19; ++angle)
            table.counts[energy][angle] = counts[particleMode][energy][angle];
    currentPlotData.angularTables.push_back(table);
    return int(currentPlotData.angularTables.size()) - 1;
}

int paceRegisterGammaEnergyTable(const QString &label, const int totalGammaCounts[33])
{
    PaceAngularTablePlotData table;
    table.label = label;
    // PACE gamma bins use the same 1 MeV boundaries as the printed table.
    table.energyLow = 1.0;
    table.energyStep = 1.0;
    table.particleMode = 4;
    table.energyOnly = true;
    table.fitKind = PaceFitNone;
    for (int energy = 0; energy < 33; ++energy)
        table.counts[energy][1] = totalGammaCounts[energy];
    currentPlotData.angularTables.push_back(table);
    return int(currentPlotData.angularTables.size()) - 1;
}

void paceRegisterCMSpectra(const int counts[5][52], const double averageEnergy[5])
{
    currentPlotData.cmSpectra = PaceCMSpectraPlotData();
    for (int particle = 0; particle < 5; ++particle)
    {
        currentPlotData.cmSpectra.averageEnergy[particle] = averageEnergy[particle];
        for (int bin = 0; bin < 52; ++bin)
            currentPlotData.cmSpectra.counts[particle][bin] = counts[particle][bin];
    }
    currentPlotData.cmSpectra.valid = true;
}

PacePlotData pacePlotDataSnapshot()
{
    return currentPlotData;
}

void openPaceAngularPlotWindow(QWidget *parent,
                               const PacePlotData &data,
                               bool plotAllTables,
                               int tableIndex,
                               PaceAngularPlotKind plotKind)
{
    if (data.angularTables.empty()) return;
    if (!plotAllTables && (tableIndex < 0 || tableIndex >= int(data.angularTables.size()))) return;

    QDialog *dialog = new QDialog(parent);
    dialog->setAttribute(Qt::WA_DeleteOnClose, true);
    dialog->setWindowIcon(QIcon(":/pace4.png"));
    dialog->setWindowTitle(plotAllTables
                               ? "PACE4: Angular distribution plots (all) - " + plotKindLabel(plotKind)
                               : "PACE4: Angular distribution plot - " + plotKindLabel(plotKind));
    dialog->setModal(false);
    QVBoxLayout *mainLayout = new QVBoxLayout(dialog);
    QHBoxLayout *topRow = new QHBoxLayout;
    QLabel *heading = new QLabel("<h2 style='color:#1d4f91; margin:0;'>Angular distribution plots: " +
                                 plotKindLabel(plotKind).toHtmlEscaped() + "</h2>");
    topRow->addWidget(heading);
    topRow->addStretch();
    mainLayout->addLayout(topRow);

    QScrollArea *scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    QWidget *container = new QWidget;
    QVBoxLayout *plots = new QVBoxLayout(container);
    int shown = 0;
    for (int i = 0; i < int(data.angularTables.size()); ++i)
    {
        if (!plotAllTables && i != tableIndex) continue;
        if (data.angularTables[i].energyOnly && plotKind != PacePlotCountsVsEnergy) continue;
        QLabel *label = new QLabel("<b>" + data.angularTables[i].label.toHtmlEscaped() + "</b>");
        label->setWordWrap(true);
        plots->addWidget(label);
        QWidget *plot = plotKind == PacePlotEnergyVsTheta
                            ? static_cast<QWidget *>(new PaceEnergyThetaPlot(data.angularTables[i]))
                            : static_cast<QWidget *>(new PaceOneDimensionalPlot(data.angularTables[i], plotKind));
        plots->addWidget(plot);
        shown++;
    }
    if (shown == 0)
    {
        delete dialog;
        return;
    }
    plots->addStretch();
    scroll->setWidget(container);
    mainLayout->addWidget(scroll);
    addSaveButton(topRow, dialog, container);
    dialog->resize(1080, std::clamp(220 + shown * 620, 700, 920));
    dialog->show();
}

void openPaceCMSpectraPlotWindow(QWidget *parent, const PacePlotData &data)
{
    if (!data.cmSpectra.valid) return;
    QDialog *dialog = new QDialog(parent);
    dialog->setAttribute(Qt::WA_DeleteOnClose, true);
    dialog->setWindowIcon(QIcon(":/pace4.png"));
    dialog->setWindowTitle("PACE4: C.M. spectra histograms");
    dialog->setModal(false);
    QVBoxLayout *mainLayout = new QVBoxLayout(dialog);
    QHBoxLayout *topRow = new QHBoxLayout;
    QLabel *heading = new QLabel("<h2 style='color:#1d4f91; margin:0;'>C.M. spectra histograms</h2>"
                                 "<p style='margin:2px 0 0 0;'><em>one histogram per emitted particle</em></p>");
    topRow->addWidget(heading);
    topRow->addStretch();
    mainLayout->addLayout(topRow);
    PaceCMSpectraPlot *plot = new PaceCMSpectraPlot(data.cmSpectra);
    mainLayout->addWidget(plot);
    addSaveButton(topRow, dialog, plot);
    dialog->resize(1200, 880);
    dialog->show();
}
