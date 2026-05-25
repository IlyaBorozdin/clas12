#pragma once

#include "iHistogramAnalyzer.h"

class BaseHistogramAnalyzer : public IHistogramAnalyzer {
public:
    BaseHistogramAnalyzer() = default;
    virtual ~BaseHistogramAnalyzer() = default;

    double getIntegral(TH1F* h, double xmin, double xmax) const override {
        int bin_min = h->FindBin(xmin);
        int bin_max = h->FindBin(1.05);

        double integral = 0.0;

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            integral += h->GetBinContent(bin);
        }

        integral /= h->GetBinWidth(1);
    }

    int countEntries(TH1F* h, double xmin, double xmax) const override {
        if (xmin > xmax) return 0;

        int bin_min = h->FindBin(xmin);
        int bin_max = h->FindBin(xmax);

        int events = 0;
        for (int bin = bin_min; bin <= bin_max; ++bin) {
            events += h->GetBinContent(bin);
        }

        return events;
    }

    void findPeak(TH1F* h, double xmin, double xmax, double& pos, double& val) const {
        int bin_min = h->FindBin(xmin);
        int bin_max = h->FindBin(1.05);

        int maxBin = bin_min;
        val = h->GetBinContent(bin_min);

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            double value = h->GetBinContent(bin);
            if (value > val) {
                val = value;
                maxBin = bin;
            }
        }

        pos = h->GetBinCenter(maxBin);
    }

    void getGlobalMax(TH1F* hist, double xmin, double xmax, double& pos, double& val) const {
        int maxBin = hist->GetMaximumBin();
        val = hist->GetBinContent(maxBin);
        pos = hist->GetBinCenter(maxBin);
    }

    void constrainToAxis(TH1F* hist, double& low, double& high) const override {
        if (!hist) return;

        double x_min = hist->GetXaxis()->GetXmin();
        double x_max = hist->GetXaxis()->GetXmax();

        if (low < x_min) low = x_min;
        if (high > x_max) high = x_max;

        // Дополнительная проверка на здравый смысл:
        // Если после коррекции границы перепутались, можно либо бросить ошибку,
        // либо схлопнуть их в одну точку.
        if (low > high) low = high; 
    }
};