#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

class BinOptimizer {
public:
    BinOptimizer(const std::vector<double>& data, double mu = 1.0, double alpha = 1.0)
        : rawData(data), nEntries(data.size()),
          paramMu(mu), paramAlpha(alpha) {
        sortedData = rawData;
        std::sort(sortedData.begin(), sortedData.end());
        estimateRange();
    }

    int findOptimalBinCount() {
        auto lower = std::lower_bound(sortedData.begin(), sortedData.end(), xMin);
        auto upper = std::upper_bound(sortedData.begin(), sortedData.end(), xMax);
        std::vector<double> filtered(lower, upper);

        int bestBins = getMinBins();
        double bestScore = std::numeric_limits<double>::infinity();

        for (int bins = getMinBins(); bins <= getMaxBins(); ++bins) {
            double score = evaluateObjective(filtered, bins);
            if (score < bestScore) {
                bestScore = score;
                bestBins = bins;
            }
        }

        return bestBins;
    }

    double getXMin() const { return xMin; }
    double getXMax() const { return xMax; }

private:
    const std::vector<double>& rawData;
    std::vector<double> sortedData;
    size_t nEntries;

    double xMin;
    double xMax;

    double paramMu;
    double paramAlpha;

    static constexpr double ABS_MIN = 0.80;
    static constexpr double XMIN_CLIP = 0.70;
    static constexpr double ABS_MAX = 1.10;
    static constexpr int ENTRIES_MIN = 3;

    int getMinBins() const { return 5; }
    int getMaxBins() const { return std::max(static_cast<int>(2 * std::sqrt(2 * nEntries)), getMinBins()); }

    void estimateRange() {
        if (nEntries < ENTRIES_MIN) {
            xMin = ABS_MIN;
            xMax = ABS_MAX;

            return;
        }

        size_t lowIndex  = static_cast<size_t>(std::ceil(0.050 * nEntries));
        size_t highIndex = static_cast<size_t>(std::floor(0.950 * nEntries));

        xMin = std::max(std::min(ABS_MIN, sortedData[lowIndex]), XMIN_CLIP);
        xMax = std::max(ABS_MAX, sortedData[highIndex]);
    }

    double evaluateObjective(const std::vector<double>& data, int binCount) {
        std::vector<int> counts(binCount, 0);
        double binWidth = (xMax - xMin) / binCount;

        size_t i = 0;
        double binLeft = xMin;
        double binRight = xMin + binWidth;
        int binIndex = 0;

        while (i < data.size() && binIndex < binCount) {
            int count = 0;
            while (i < data.size() && data[i] < binRight) {
                ++count;
                ++i;
            }
            counts[binIndex] = count;
            ++binIndex;
            binLeft = binRight;
            binRight += binWidth;
        }

        double contrast = computeContrast(counts);
        double penalty = binPenalty(binCount);
        return contrast + penalty;
    }

    double computeContrast(const std::vector<int>& h) {
        double totalContrast = 0.0;

        for (size_t i = 1; i + 1 < h.size(); ++i) {
            if (h[i] >= h[i - 1] && h[i] >= h[i + 1]) {
                double meanNeighbor = 0.5 * (h[i - 1] + h[i + 1]);
                if (h[i] > 0) {
                    double ci = (h[i] - meanNeighbor) / h[i];
                    totalContrast += ci;
                }
            }
        }

        if (h[0] >= h[1] && h[0] > 0) {
            double ci = (h[0] - 0.5 * h[1]) / h[0];
            totalContrast += ci;
        }

        if (h[h.size() - 1] >= h[h.size() - 2] && h[h.size() - 1] > 0) {
            double ci = (h[h.size() - 1] - 0.5 * h[h.size() - 2]) / h[h.size() - 1];
            totalContrast += ci;
        }

        return totalContrast;
    }

    double binPenalty(int binCount) {
        return paramMu / std::pow(binCount, paramAlpha);
    }
};
