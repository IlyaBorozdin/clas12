#pragma once

#include "fittingStrategyFactory.h"

class HistogramFitterFlexible : public FittingStrategy {
public:
    TF1* fit(TH1F* hist, int i, int j, int k, int l) const override { 
        FitType type = chooseFitType(i, j, k, l);
        std::unique_ptr<FittingStrategy> strategy = FittingStrategyFactory::create(type);
        return strategy->fit(hist, i, j, k, l); 
    }
    
    double getDownEdge(TH1F* hist, int i, int j, int k, int l) const override {
        FitType type = chooseFitType(i, j, k, l);
        std::unique_ptr<FittingStrategy> strategy = FittingStrategyFactory::create(type);
        return strategy->getDownEdge(hist, i, j, k, l);
    }

    double getUpEdge(TH1F* hist, int i, int j, int k, int l) const override {
        FitType type = chooseFitType(i, j, k, l);
        std::unique_ptr<FittingStrategy> strategy = FittingStrategyFactory::create(type);
        return strategy->getUpEdge(hist, i, j, k, l);
    }

    double getDeltaPeak(TH1F* hist, int i, int j, int k, int l) const override {
        FitType type = chooseFitType(i, j, k, l);
        std::unique_ptr<FittingStrategy> strategy = FittingStrategyFactory::create(type);
        return strategy->getDeltaPeak(hist, i, j, k, l);
    }

private:
    // static constexpr double W_RESOLUTION_THRESHOLD = 1.35;
    // static constexpr double PHI_RESOLUTION_THRESHOLD_DOWN = 0.3 * M_PI;
    // static constexpr double PHI_RESOLUTION_THRESHOLD_UP = 0.7 * M_PI;

    /*
    double getW(int j) const {
        return W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W;
    }

    double getPhi(int l) const {
        return (l + 0.5) * 2 * M_PI / NUMBER_PHI;
    }

    FitType chooseFitType(int i, int j, int k, int l) const {
        // Примерная логика выбора (можно улучшить/расширить):
        if (getW(j) > W_RESOLUTION_THRESHOLD &&
            getPhi(l) > PHI_RESOLUTION_THRESHOLD_DOWN &&
            getPhi(l) < PHI_RESOLUTION_THRESHOLD_UP) return FitType::TwoPeaks;
        return FitType::SinglePeak;
    }
    */

    bool isSiglePeakMode(int i, int j, int k, int l) const {
        convertVariableToNormal(i, j, k, l);

        if (j >= 1 && j <= 5) return true;
        if (j >= 6 && j <= 10 && l >= 4 && l <= 7) return true;
        if (j >= 11 && j <= 12 && l >= 5 && l <= 6) return true;
        return false;
    }

    bool isTwoPeakMode(int i, int j, int k, int l) const {
        convertVariableToNormal(i, j, k, l);

        if (j >= 6 && j <= 10 && (l <= 3 || l >= 8)) return true;
        if (j >= 11 && j <= 12 && (l <= 4 || l >= 7)) return true;
        if (j >= 13 && j <= 22) return true;
        return false;
    }

    FitType chooseFitType(int i, int j, int k, int l) const {
        if (isSiglePeakMode(i, j, k, l)) return FitType::SinglePeak;
        if (isTwoPeakMode(i, j, k, l)) return FitType::TwoPeaks;
        throw std::runtime_error("Unknown fitting strategy");
    }
};
