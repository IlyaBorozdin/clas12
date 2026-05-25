#pragma once

#include "fitterStrategy.h"

class BaseFitter : public FitterStrategy {
public:
    using FitterStrategy::FitterStrategy;

    virtual double getDownEdge(TH1F* hist, int i, int j, int k, int l) const = 0;
    virtual double getUpEdge(TH1F* hist, int i, int j, int k, int l) const = 0;
    virtual double getDeltaPeak(TH1F* hist, int i, int j, int k, int l) const = 0;
    
    void getNeutronPeak(TH1F* hist, int i, int j, int k, int l, double& maxPosition, double& maxValue) const {
        double downEdge = getDownEdge(hist, i, j, k, l);

        int bin_min = hist->FindBin(downEdge);
        int bin_max = hist->FindBin(1.05);

        int maxBin = bin_min;
        maxValue = hist->GetBinContent(bin_min);

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            double value = hist->GetBinContent(bin);
            if (value > maxValue) {
                maxValue = value;
                maxBin = bin;
            }
        }

        maxPosition = hist->GetBinCenter(maxBin);
    }

    void getNeutronIntegral(TH1F* hist, int i, int j, int k, int l, double& integral) const {
        double downEdge = getDownEdge(hist, i, j, k, l);

        int bin_min = hist->FindBin(downEdge);
        int bin_max = hist->FindBin(1.05);

        integral = 0.0;

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            integral += hist->GetBinContent(bin);
        }

        integral /= hist->GetBinWidth(1);
    }

    void getMaxBin(TH1F* hist, double& maxPosition, double& maxValue) const {
        int maxBin = hist->GetMaximumBin();
        maxValue = hist->GetBinContent(maxBin);
        maxPosition = hist->GetBinCenter(maxBin);
    }

    int getFitRegionEvents(TH1F* hist, int i, int j, int k, int l) const {
        double downEdge = getDownEdge(hist, i, j, k, l);
        double upEdge = getUpEdge(hist, i, j, k, l);

        if (downEdge > upEdge) return 0;

        int bin_min = hist->FindBin(downEdge);
        int bin_max = hist->FindBin(upEdge);

        int events = 0;
        for (int bin = bin_min; bin <= bin_max; ++bin) {
            events += hist->GetBinContent(bin);
        }

        return events;
    }
};