#pragma once

#include "singlePeakFittingStrategy.h"

class TwoPeakFittingStrategy : public FittingStrategy {
public:
    TF1* fit(TH1F* hist, int i, int j, int k, int l) const override {
        /**
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + [4]*(x - [5])*(x - [6]) + [7]*exp(-0.5*((x-[8])/[9])^2)";

        double upEdge = getUpEdge(hist, i, j, k, l);
        double downEdge = getDownEdge(hist, i, j, k, l);
        double deltaPeak = getDeltaPeak(hist, i, j, k, l);

        TF1* fitFunc = new TF1("Two Peak Fit Function", strFitFunc.c_str(), downEdge, upEdge);

        SinglePeakFittingStrategy primaryFitter;
        TF1* primaryFitFunc = primaryFitter.fit(hist, i, j, k, l);

        if (!primaryFitFunc) return nullptr;

        for (int f = 0; f < 7; f++) {
            fitFunc->SetParameter(f, primaryFitFunc->GetParameter(f));
        }

        fitFunc->SetParameter(7, 0.8 * hist->GetBinContent(hist->FindBin(deltaPeak)));
        fitFunc->FixParameter(8, deltaPeak);
        fitFunc->SetParameter(9, 0.03);

        fitFunc->SetParLimits(0, 0.7 * primaryFitFunc->GetParameter(0), 1.3 * primaryFitFunc->GetParameter(0));
        fitFunc->SetParLimits(1, 0.9 * primaryFitFunc->GetParameter(1), 1.1 * primaryFitFunc->GetParameter(1));
        fitFunc->SetParLimits(2, 0.9 * primaryFitFunc->GetParameter(2), 1.1 * primaryFitFunc->GetParameter(2));
        fitFunc->SetParLimits(3, 0.9 * primaryFitFunc->GetParameter(3), 1.1 * primaryFitFunc->GetParameter(3));
        fitFunc->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(deltaPeak)));
        fitFunc->SetParLimits(9, 0.01, 0.1);

        hist->Fit(fitFunc, "RQ");

        delete primaryFitFunc;
        return fitFunc;
        */

        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + TMath::Max(0., [4]*(x - [5])*(x - [6])) + [7]*exp(-0.5*((x-[8])/[9])^2)";

        double upEdge = getUpEdge(hist, i, j, k, l);
        double downEdge = getDownEdge(hist, i, j, k, l);
        double deltaPeak = getDeltaPeak(hist, i, j, k, l);

        double maxPosition, maxValue;
        getNeutronPeak(hist, i, j, k, l, maxPosition, maxValue);

        TF1* fitFunc = new TF1("Two Peak Fit Function", strFitFunc.c_str(), downEdge, upEdge);

        // Инициализация параметров первой функции
        fitFunc->SetParameter(0, maxValue);
        fitFunc->SetParameter(1, maxPosition);
        fitFunc->SetParameter(2, 0.03);
        fitFunc->SetParameter(3, 1.00);

        if (isEmptyBack(i, j, k, l)) {
            fitFunc->FixParameter(4, 0.00);
            fitFunc->FixParameter(5, 0.00);
            fitFunc->FixParameter(6, 0.00);
            
        } else {
            double x1 = 0.00;
            double x2 = 0.00;

            getZerosParabola(i, j, k, l, x1, x2);

            // double meanHeight = computeMeanBackground(hist, 1.00, x2);
            // double a = computeQuadraticCoeff(meanHeight, x1, x2);

            // fitFunc->SetParameter(4, a / 5.0);
            if (x2 > 0.8) fitFunc->SetParameter(4, -13.00);
            else fitFunc->SetParameter(4, 6.00);
            
            fitFunc->FixParameter(5, x1);
            fitFunc->FixParameter(6, x2);
        }

        fitFunc->SetParameter(7, 0.8 * hist->GetBinContent(hist->FindBin(deltaPeak)));
        fitFunc->FixParameter(8, deltaPeak);
        fitFunc->SetParameter(9, 0.03);

        // Ограничения для первой аппроксимации
        fitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        fitFunc->SetParLimits(1, 0.85, 1.05);
        fitFunc->SetParLimits(2, 0.0085, 0.085);
        fitFunc->SetParLimits(3, 0.33, 3.00);
        fitFunc->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(deltaPeak)));
        fitFunc->SetParLimits(9, 0.01, 0.1);
        // fitFunc->SetParLimits(5, 0.60, 0.95);
        // fitFunc->SetParLimits(6, 0.60, 4.00);

        hist->Fit(fitFunc, "RQ");

        return fitFunc;
    }

    double getDownEdge(TH1F* hist, int i, int j, int k, int l) const override {
        double downEdge = 0.8;
        double x_min = hist->GetXaxis()->GetXmin();

        if (downEdge < x_min) downEdge = x_min;
        return downEdge;
    }

    double getUpEdge(TH1F* hist, int i, int j, int k, int l) const override {
        double upEdge = MM_UP_EDGES[j * 18 / NUMBER_W];
        double x_max = hist->GetXaxis()->GetXmax();

        if (upEdge > x_max) upEdge = x_max;
        return upEdge;
    }

    double getDeltaPeak(TH1F* hist, int i, int j, int k, int l) const override {
        return DELTA_PEAK_POSITIONS[j * 18 / NUMBER_W];
    }

private:
    static constexpr double MM_UP_EDGES[18] = {
        0.99, 1.04, 1.06, 1.09, 1.09, 1.21, 1.25, 1.30, 1.34,
        1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38
    };

    static constexpr double DELTA_PEAK_POSITIONS[18] = {
        1.000, 1.000, 1.000, 1.000, 1.000, 1.088, 1.134, 1.169, 1.201, 
        1.223, 1.239, 1.230, 1.227, 1.233, 1.216, 1.220, 1.220, 1.220
    };

    bool isEmptyBack(int i, int j, int k, int l) const {
        ++i; ++j; ++k; ++l;
        if (j == 6 && (k <= 2 || k >= 8)) return true;
        if (j >= 7 && j <= 8 && (k >= 8 && k <= 10)) return true;
        if (j >= 9 && j <= 10 && (k >= 9 && k <= 10)) return true;
        if (j == 11 && k == 10) return true;
        if (j >= 12 && j <= 13 && (k >= 2 && k <= 6) && (l >= 5 && l <= 6)) return true;
        if (j == 14 && (k >= 3 && k <= 5) && (l >= 5 && l <= 6)) return true;
        if (j >= 15 && j <= 18 && k == 2 && (l >= 5 && l <= 6)) return true;
        return false;
    }

    void getZerosParabola(int i, int j, int k, int l, double& x1, double& x2) const {
        ++i; ++j; ++k; ++l;
        if (j == 6 && k >= 3 && k <= 4 && (l <= 3 || l >= 8)) { x1 = 0.80; x2 = 1.20; }
        else if (j == 6 && k >= 3 && k <= 4 && (l >= 4 && l <= 7)) { x1 = 0.78; x2 = 1.12; }
        else if (j == 6 && k >= 5 && k <= 6 && (l <= 3 || l >= 8)) { x1 = 0.82; x2 = 1.18; }
        else if (j == 6 && k >= 5 && k <= 6 && (l >= 4 && l <= 7)) { x1 = 0.78; x2 = 1.12; }
        else if (j == 6 && k == 7 && (l <= 3 || l >= 8)) { x1 = 0.84; x2 = 1.16; }
        else if (j == 6 && k == 7 && (l >= 4 && l <= 7)) { x1 = 0.80; x2 = 1.10; }

        else if (j == 7 && k >= 1 && k <= 7 && (l <= 3 || l >= 8)) { x1 = 0.85; x2 = 1.27; }
        else if (j == 7 && k >= 1 && k <= 7 && (l >= 4 && l <= 7)) { x1 = 0.85; x2 = 1.15; }

        else if (j == 8 && k >= 1 && k <= 7 && (l <= 3 || l >= 8)) { x1 = 0.84; x2 = 1.30; }
        else if (j == 8 && k >= 1 && k <= 7 && (l >= 4 && l <= 7)) { x1 = 0.86; x2 = 1.14; }

        else if (j == 9 && k >= 1 && k <= 8 && (l <= 3 || l >= 8)) { x1 = 0.83; x2 = 1.33; }
        else if (j == 9 && k >= 1 && k <= 8 && (l >= 4 && l <= 7)) { x1 = 0.86; x2 = 1.14; }

        else if (j == 10 && k >= 1 && k <= 8 && (l <= 4 || l >= 7)) { x1 = 0.83; x2 = 1.37; }
        else if (j == 10 && k >= 1 && k <= 8 && (l >= 5 && l <= 6)) { x1 = 0.86; x2 = 1.14; }

        else if (j == 11 && k >= 1 && k <= 9 && (l <= 4 || l >= 7)) { x1 = 0.83; x2 = 1.37; }
        else if (j == 11 && k >= 1 && k <= 9 && (l >= 5 && l <= 6)) { x1 = 0.86; x2 = 1.14; }

        else if (j == 12 && k == 1 ||
                 j == 12 && k == 2 && (l >= 2 && l <= 4 || l >= 7 && l <= 9) ||
                 j == 12 && k >= 3 && k <= 4 && (l >= 3 && l <= 4 || l >= 7 && l <= 8) ||
                 j == 12 && k >= 5 && k <= 6 && (l >= 3 && l <= 4 || l >= 7 && l <= 8) ||
                 j == 12 && k >= 7 && k <= 8 && (l >= 3 && l <= 7) ||
                 j == 12 && k >= 9 && k <= 10 ||
                 
                 j == 13 && k == 1 ||
                 j == 13 && k == 2 && (l >= 2 && l <= 4 || l >= 7 && l <= 9) ||
                 j == 13 && k == 3 && (l >= 3 && l <= 4 || l >= 7 && l <= 8) ||
                 j == 13 && k >= 4 && k <= 6 && (l == 4 || l == 7) ||
                 j == 13 && k >= 7 && k <= 8 && (l >= 5 && l <= 6) ||
                 
                 j == 14 && k == 1 ||
                 j == 14 && k == 2 && (l >= 2 && l <= 9) ||
                 j == 14 && k == 3 && (l >= 3 && l <= 4 || l >= 7 && l <= 8) ||
                 j == 14 && k >= 4 && k <= 5 && (l == 4 && l == 7) ||
                 j == 14 && k >= 6 && k <= 7 && (l >= 5 && l <= 6) ||
                 j == 14 && k == 10 ||
                 
                 j >= 15 && j <= 18 && k == 1 ||
                 j >= 15 && j <= 17 && k == 2 && (l >= 2 && l <= 4 || l >= 7 && l <= 9) ||
                 j >= 15 && j <= 18 && k >= 3 && k <= 5 && (l >= 4 && l <= 7) ||
                 j >= 15 && j <= 17 && k == 6 && (l >= 5 && l <= 6) ||
                 j >= 15 && j <= 18 && k == 10 ||
                 
                 j == 18 && k == 2 && (l >= 3 && l <= 4 || l == 7)) { x1 = 0.83; x2 = 1.37; }

        else if (j == 13 && k >= 9 && k <= 10) { x1 = 0.80; x2 = 1.40; }

        else if (j == 12 && k == 2 && (l <= 1 || l >= 10) ||
                 j == 12 && k >= 3 && k <= 4 && (l <= 2 || l >= 9) ||
                 j == 12 && k >= 5 && k <= 6 && (l <= 2 || l >= 9) ||
                 j == 12 && k >= 7 && k <= 8 && (l <= 2 || l >= 8) ||
                 
                 j == 13 && k == 2 && (l <= 1 || l >= 10) ||
                 j == 13 && k == 3 && (l <= 2 || l >= 9) ||
                 j == 13 && k >= 4 && k <= 6 && (l <= 3 || l >= 8) ||
                 j == 13 && k >= 7 && k <= 8 && (l <= 4 || l >= 7) ||
                 
                 j == 14 && k == 2 && (l <= 1 || l >= 10) ||
                 j == 14 && k == 3 && (l <= 2 || l >= 9) ||
                 j == 14 && k >= 4 && k <= 5 && (l <= 3 || l >= 8) ||
                 j == 14 && k >= 6 && k <= 7 && (l <= 4 || l >= 7) ||
                 j == 14  && k >= 8 && k <= 9 ||
                 
                 j >= 15 && j <= 17 && k == 2 && (l <= 1 || l >= 10) ||
                 j >= 15 && j <= 18 && k >= 3 && k <= 5 && (l <= 3 || l >= 8) ||
                 j >= 15 && j <= 17 && k == 6 && (l <= 4 || l >= 7) ||
                 j >= 15 && j <= 17 && k >= 7 && k <= 9 ||
                 
                 j == 18 && k == 2 && (l <= 2 || l >= 8) ||
                 j == 18 && k == 6) { x1 = 0.70; x2 = 0.70; }
    }
};
