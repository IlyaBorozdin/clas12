#pragma once

#include <TH1F.h>
#include <TF1.h>
#include "iParameterProvider.h"
#include "iHistogramAnalyzer.h"

class FitterStrategy {
protected:
    IParameterProvider* m_provider; // Провайдер, задающий "числа"
    IHistogramAnalyzer* m_analyzer;

public:
    FitterStrategy(IParameterProvider* provider, IHistogramAnalyzer* analyzer) 
        : m_provider(provider), m_analyzer(analyzer) {}
    virtual ~FitterStrategy() = default;

    // Главный метод, который должен реализовать каждый конкретный фиттер
    virtual TF1* fit(TH1F* hist, int i, int j, int k, int l) const = 0;

    virtual double calculateYield(TF1* func, double binWidth) const = 0;
    virtual double calculateYieldError(TF1* func, double binWidth) const = 0;

    // const IParameterProvider* getProvider() const {
    //     return m_provider;
    // }

    const IHistogramAnalyzer* getAnalyzer() const {
        return m_analyzer;
    }

    void getBounds(TH1F* h, int i, int j, int k, int l, double& low, double& high) const {
        m_provider->getBounds(i, j, k, l, low, high);
        m_analyzer->constrainToAxis(h, low, high);
    }
};