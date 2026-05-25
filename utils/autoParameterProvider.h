#pragma once

#include "iParameterProvider.h"
#include "autoHistogramAnalyzer.h"
#include "../MM_project_const.h"

class AutoParameterProvider : public IParameterProvider {
public:
    AutoParameterProvider(AutoHistogramAnalyzer* analyzer)
        : IParameterProvider(analyzer) {}

    FitParameter getParameter(const std::string& name, int i, int j, int k, int l, TH1F* hist = nullptr) const override {
        if (!hist) return FitParameter();

        
    }

    bool hasParameter(const std::string& name, int i, int j, int k, int l) const override {

    }

    void getBounds(int i, int j, int k, int l, double& low, double& high) const override {
        low = 0.65;
        high = 1.80;
    }
};