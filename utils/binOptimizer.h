#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

struct WeightedValue {
    double value;
    double weight;

    WeightedValue(double v = 0.0, double w = 1.0)
        : value(v), weight(w) {}
};

class BinOptimizer {
public:
    // Старый интерфейс: массив значений без весов
    BinOptimizer(const std::vector<double>& data, double mu = 1.0, double alpha = 1.0)
        : nEntries(data.size()), paramMu(mu), paramAlpha(alpha) {
        rawData.reserve(data.size());
        for (double v : data) {
            rawData.emplace_back(v, 1.0); // все веса = 1
        }
        sortedData = rawData;
        sortData();
        estimateRange();
    }

    // Новый интерфейс: массив WeightedValue
    BinOptimizer(const std::vector<WeightedValue>& data, double mu = 1.0, double alpha = 1.0)
        : rawData(data), sortedData(data), nEntries(data.size()), paramMu(mu), paramAlpha(alpha) {
        sortData();
        estimateRange();
    }

    int findOptimalBinCount() {
        std::vector<WeightedValue> filtered;

        std::copy_if(
            sortedData.begin(), sortedData.end(),
            std::back_inserter(filtered),
            [this](const WeightedValue& wv) {
                return wv.value >= this->xMin && wv.value <= this->xMax;
            }
        );

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
    std::vector<WeightedValue> rawData;
    std::vector<WeightedValue> sortedData;
    size_t nEntries;

    double xMin;
    double xMax;

    double paramMu;
    double paramAlpha;

    static constexpr double ABS_MIN = 0.80;
    static constexpr double XMIN_CLIP = 0.65;

    static constexpr double ABS_MAX = 1.10;
    static constexpr double XMAX_CLIP = 1.80;

    int getMinBins() const { return 5; }
    int getMaxBins() const { return std::max(static_cast<int>(2 * std::sqrt(2 * nEntries)), getMinBins()); }

    void sortData() {
        std::sort(sortedData.begin(), sortedData.end(),
                  [](const WeightedValue& a, const WeightedValue& b) {
                      return a.value < b.value;
                  });
    }

    void estimateRange() {
        if (nEntries == 0) {
            xMin = ABS_MIN;
            xMax = ABS_MAX;
            return;
        }

        xMin = std::max(std::min(ABS_MIN, sortedData.front().value), XMIN_CLIP);
        xMax = std::max(ABS_MAX, std::min(XMAX_CLIP, sortedData.back().value));
    }

    double evaluateObjective(const std::vector<WeightedValue>& data, int binCount) {
        std::vector<double> binWeights(binCount, 0.0);
        double binWidth = (xMax - xMin) / binCount;

        size_t i = 0;
        int binIndex = 0;
        double binRight = xMin + binWidth;

        while (i < data.size() && binIndex < binCount) {
            while (i < data.size() && data[i].value < binRight) {
                binWeights[binIndex] += data[i].weight; // добавляем вес
                ++i;
            }
            ++binIndex;
            binRight += binWidth;
        }

        double contrast = computeContrast(binWeights);
        double penalty = binPenalty(binCount);
        return contrast + penalty;
    }

    double computeContrast(const std::vector<double>& h) {
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

        if (h.size() >= 2) {
            if (h[0] >= h[1] && h[0] > 0) {
                double ci = (h[0] - 0.5 * h[1]) / h[0];
                totalContrast += ci;
            }
            if (h[h.size() - 1] >= h[h.size() - 2] && h[h.size() - 1] > 0) {
                double ci = (h[h.size() - 1] - 0.5 * h[h.size() - 2]) / h[h.size() - 1];
                totalContrast += ci;
            }
        }

        return totalContrast;
    }

    double binPenalty(int binCount) {
        return paramMu / std::pow(binCount, paramAlpha);
    }
};
