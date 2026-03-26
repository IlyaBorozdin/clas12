#pragma once

#include "../MM_project_const.h"

#include "TH1F.h"
#include "TF1.h"

class FittingStrategy {
public:
    virtual ~FittingStrategy() = default;

    virtual TF1* fit(TH1F* hist, int i, int j, int k, int l) const = 0;

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

    double computeMeanBackground(TH1F* hist, double x1, double x2) const {
        double x_min = hist->GetXaxis()->GetXmin();
        double x_max = hist->GetXaxis()->GetXmax();

        if (x1 < x_min) x1 = x_min;
        if (x2 > x_max) x2 = x_max;
        if (x1 > x2) return 0.0; 

        int bin_min = hist->FindBin(x1);
        int bin_max = hist->FindBin(x2);

        double sum = 0.0;
        int count = bin_max - bin_min + 1;
        for (int bin = bin_min; bin <= bin_max; ++bin) {
            sum += hist->GetBinContent(bin);
        }

        return sum / (bin_max - bin_min + 1);
    }

    double computeQuadraticCoeff(double meanHeight, double x1, double x2) const {
        double dx = x2 - x1;
        if (dx <= 0.0) return 0.0;
        return -4.0 * meanHeight / (dx * dx);
    }

    void convertNumbersToNormalOrder(int& i, int& j, int& k, int& l) const {
        ++i; ++j; ++k; ++l;
    }

    void convertNumbersToVariable(int& i, int& j, int& k, int& l, double& q2, double& w, double& cos_theta, double& phi) const {
        q2 = (STEPS_Q2[2 * i] + STEPS_Q2[2 * i + 1]) / 2;
        w = (0.5 + j) * STEP_W + W_MIN;
        cos_theta = (0.5 + k) * STEP_COS_THETA - 1;
        phi = (0.5 + l) * STEP_PHI;
    }
};
