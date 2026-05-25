#pragma once

#include "baseHistogramAnalyzer.h"

class AutoHistogramAnalyzer : public BaseHistogramAnalyzer {
public:
    AutoHistogramAnalyzer() = default;
    virtual ~AutoHistogramAnalyzer() = default;

    double getIntegral(TH1F* h, double xmin, double xmax) const override {

    }


    void findPeak(TH1F* h, double xmin, double xmax, double& pos, double& val) const {

    }
};