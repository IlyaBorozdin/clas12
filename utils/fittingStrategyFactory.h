#pragma once

#include <memory>

#include "twoPeakFittingStrategy.h"

enum class FitType {
    SinglePeak,
    TwoPeaks,
    EmptyBack,
};

class FittingStrategyFactory {
public:
    static std::unique_ptr<FittingStrategy> create(FitType type) {
        switch (type) {
            case FitType::SinglePeak: return std::make_unique<SinglePeakFittingStrategy>();
            case FitType::TwoPeaks: return std::make_unique<TwoPeakFittingStrategy>();
            default: throw std::runtime_error("Unknown fitting strategy");
        }
    }
};
