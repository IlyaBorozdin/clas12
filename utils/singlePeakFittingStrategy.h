#pragma once

#include "fittingStrategy.h"
#include "../MM_project_const.h"

class SinglePeakFittingStrategy : public FittingStrategy {
public:
    TF1* fit(TH1F* hist, int i, int j, int k, int l) const override {
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + TMath::Max(0., [4]*(x - [5])*(x - [6])) + [7]*exp(-0.5*((x-[8])/[9])^2)";

        double upEdge = getUpEdge(hist, i, j, k, l);
        double downEdge = getDownEdge(hist, i, j, k, l);

        double maxPosition, maxValue;
        getNeutronPeak(hist, i, j, k, l, maxPosition, maxValue);

        TF1* fitFunc = new TF1("Single Peak Fit Function", strFitFunc.c_str(), downEdge, upEdge);

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
            
            fitFunc->SetParameter(4, 0.00);
            fitFunc->FixParameter(5, x1);
            fitFunc->FixParameter(6, x2);
        }

        fitFunc->FixParameter(7, 0.00);
        fitFunc->FixParameter(8, 0.00);
        fitFunc->FixParameter(9, 0.03);

        // Ограничения для первой аппроксимации
        fitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        fitFunc->SetParLimits(1, downEdge, 1.05);
        fitFunc->SetParLimits(2, 0.0085, 0.085);
        fitFunc->SetParLimits(3, 0.33, 3.00);
        // fitFunc->SetParLimits(5, 0.60, 0.95);
        // fitFunc->SetParLimits(6, 0.60, 4.00);

        hist->Fit(fitFunc, "RQ");

        return fitFunc;
    }

    double getDownEdge(TH1F* hist, int i, int j, int k, int l) const override {
        double downEdge = 0.80;

        convertVariableToNormal(i, j, k, l);

        if ((j == 1 && k == 7 && (l <= 3 || l >= 8)) ||
            (j == 4 && k == 8 && l == 6)) { downEdge = 0.75; }
        else if ((j == 1 && (k >= 6 && k <= 7) && (l >= 4 && l <= 7)) ||
                 (j == 2 && k == 8) ||
                 (j == 3 && k == 8 && (l >= 4 && l <= 7)) ||
                 (j == 4 && k == 8 && l == 5) ||
                 (j == 4 && k == 9)) { downEdge = 0.70; }
        else if ((j == 1 && (k >= 8 && k <= 10)) ||
                 ((j >= 2 && j <= 3) && (k >= 8 && k <= 10)) || // В тетради при k >= 8 !!!
                 ((j >= 4 && j <= 5) && k == 10)) { downEdge = 0.65; }

        double x_min = hist->GetXaxis()->GetXmin();

        if (downEdge < x_min) downEdge = x_min;
        return downEdge;
    }

    double getUpEdge(TH1F* hist, int i, int j, int k, int l) const override {
        double upEdge = MM_UP_SHORT_EDGES[j];
        double x_max = hist->GetXaxis()->GetXmax();

        if (upEdge > x_max) upEdge = x_max;
        return upEdge;
    }

    double getDeltaPeak(TH1F* hist, int i, int j, int k, int l) const override {
        return 0.0;
    }

private:
    static constexpr double MM_UP_SHORT_EDGES[22] = {
        0.99, 1.04, 1.06, 1.07, 1.07, 1.07, 1.07, 1.07, 1.08,
        1.08, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10
    };

    bool isEmptyBack(int i, int j, int k, int l) const {
        convertVariableToNormal(i, j, k, l);

        if (j >= 1 && j <= 3) return true;
        if (j >= 4 && j <= 5 && k >= 1 && k <= 7 && l >= 4 && l <= 7) return true;
        if (j == 4 && j <= 5 && k >= 8 && k <= 10) return true;
        return false;
    }

    void getZerosParabola(int i, int j, int k, int l, double& x1, double& x2) const {
        convertVariableToNormal(i, j, k, l);

        if (j == 4 && k >= 1 && k <= 5 && (l <= 3 || l >= 8)) { x1 = 0.86; x2 = 1.10; }
        else if (j == 4 && k == 6 && (l <= 3 || l >= 8)) { x1 = 0.88; x2 = 1.08; }
        else if (j == 4 && k == 7 && (l <= 3 || l >= 8)) { x1 = 0.88; x2 = 1.06; }
        else if (j == 5 && k >= 1 && k <= 5 && (l <= 2 || l >= 9)) { x1 = 0.86; x2 = 1.16; }
        else if (j == 5 && k >= 1 && k <= 5 && (l == 3 || l == 8)) { x1 = 0.86; x2 = 1.14; }
        else if (j == 5 && k == 6 && (l <= 2 || l >= 9)) { x1 = 0.88; x2 = 1.14; }
        else if (j == 5 && k == 6 && (l == 3 || l == 8)) { x1 = 0.88; x2 = 1.12; }
        else if (j == 5 && k == 7 && (l <= 3 || l >= 8)) { x1 = 0.88; x2 = 1.10; }
    }
};