#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "../MM_project_const.h"

class HistogramFitter {
public:
    TF1* fit(TH1F* hist, int i, int j, int k, int l) {
        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
        std::string strPeakFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
        std::string strFitFunc = strPeakFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        double upEdge = getUpEdge(i, j, k, l);
        double upShortEdge = getUpShortEdge(i, j, k, l);
        double downEdge = getDownEdge(i, j, k, l);
        double downShortEdge = getDownShortEdge(i, j, k, l);
        double deltaPeak = getDeltaPeak(i, j, k, l);

        double maxPosition, maxValue;
        getNeutronPeak(hist, maxPosition, maxValue);

        TF1* fitFunc = new TF1("Fit Function", strFitFunc.c_str(), downEdge, upEdge);
        TF1* primaryFitFunc = new TF1("Primary Fit Function", strPeakFunc.c_str(), downShortEdge, upShortEdge);

        primaryFitFunc->SetParameter(0, maxValue);
        primaryFitFunc->SetParameter(1, maxPosition);
        primaryFitFunc->SetParameter(2, 0.03);
        primaryFitFunc->SetParameter(3, 1.00);
        primaryFitFunc->SetParameter(4, 0.00);
        primaryFitFunc->SetParameter(5, 0.80);
        primaryFitFunc->SetParameter(6, 1.60);

        primaryFitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        primaryFitFunc->SetParLimits(1, 0.86, 1.00);
        primaryFitFunc->SetParLimits(2, 0.0085, 0.085);
        primaryFitFunc->SetParLimits(3, 0.33, 3.00);
        primaryFitFunc->SetParLimits(5, 0.70, 0.95);
        primaryFitFunc->SetParLimits(6, 0.70, 4.00);

        hist->Fit(primaryFitFunc, "RQ");

        for (int f = 0; f < 7; f++) {
            fitFunc->SetParameter(f, primaryFitFunc->GetParameter(f));
        }
        fitFunc->SetParameter(7, 1.00);
        fitFunc->SetParameter(8, deltaPeak);
        fitFunc->SetParameter(9, 0.03);

        fitFunc->SetParLimits(0, 0.7 * primaryFitFunc->GetParameter(0), 1.3 * primaryFitFunc->GetParameter(0));
        fitFunc->SetParLimits(1, 0.9 * primaryFitFunc->GetParameter(1), 1.1 * primaryFitFunc->GetParameter(1));
        fitFunc->SetParLimits(2, 0.9 * primaryFitFunc->GetParameter(2), 1.1 * primaryFitFunc->GetParameter(2));
        fitFunc->SetParLimits(3, 0.9 * primaryFitFunc->GetParameter(3), 1.1 * primaryFitFunc->GetParameter(3));
        fitFunc->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(deltaPeak)));
        fitFunc->SetParLimits(8, 0.8 * deltaPeak, 1.2 * deltaPeak);
        fitFunc->SetParLimits(9, 0.01, 0.1);

        hist->Fit(fitFunc, "RQ");

        delete primaryFitFunc;
        return fitFunc;
    }

    double getUpEdge(int i, int j, int k, int l) const { return interpolate(MM_UP_EDGES, j); }

    double getUpShortEdge(int i, int j, int k, int l) const { return interpolate(MM_UP_SHORT_EDGES, j); }

    double getDownEdge(int i, int j, int k, int l) const { return 0.80; }

    double getDownShortEdge(int i, int j, int k, int l) const { return 0.80; }

    double getDeltaPeak(int i, int j, int k, int l) const { return interpolate(DELTA_PEAK_POSITIONS, j); }

    void getNeutronPeak(TH1F* hist, double& maxPosition, double& maxValue) {
        int bin_min = hist->FindBin(0.86);
        int bin_max = hist->FindBin(1.00);

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

private:
    static constexpr double MM_UP_EDGES[36] = {
        0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.14, 1.18, 1.23, 1.24, 1.24, 1.28, 1.30,
        1.32, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
        1.34, 1.34, 1.34, 1.34
    };

    static constexpr double MM_UP_SHORT_EDGES[36] = {
        0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.09, 1.10, 1.11, 1.12, 1.12, 1.10, 1.08, 1.08, 1.08, 1.10,
        1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10,
        1.10, 1.10, 1.10, 1.10
    };

    static constexpr double DELTA_PEAK_POSITIONS[36] = {
        0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.96, 0.96, 0.98, 1.08, 1.08, 1.12, 1.14, 1.16, 1.16, 1.18,
        1.20, 1.20, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22,
        1.22, 1.22, 1.22, 1.22
    };

    double interpolate(const double* array, int j) const {
        const int n = 36;
        double x = static_cast<double>(j) * (n - 1) / (NUMBER_W - 1);
        int i = static_cast<int>(x);
        double frac = x - i;

        if (i >= n - 1) return array[n - 1];
        return array[i] + frac * (array[i + 1] - array[i]);
    }
};

class HistogramFitterExt {
public:
    TF1* fit(TH1F* hist, int i, int j, int k, int l) {
        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
        std::string strPeakFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
        std::string strFitFunc = strPeakFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        double upEdge = getUpEdge(i, j, k, l);
        double upShortEdge = getUpShortEdge(i, j, k, l);
        double downEdge = getDownEdge(i, j, k, l);
        double downShortEdge = getDownShortEdge(i, j, k, l);
        double deltaPeak = getDeltaPeak(i, j, k, l);

        double maxPosition, maxValue;
        getNeutronPeak(hist, maxPosition, maxValue);

        TF1* fitFunc = new TF1("Fit Function", strFitFunc.c_str(), downEdge, upEdge);
        TF1* primaryFitFunc = new TF1("Primary Fit Function", strPeakFunc.c_str(), downShortEdge, upShortEdge);

        // Инициализация параметров первой функции
        primaryFitFunc->SetParameter(0, maxValue);
        primaryFitFunc->SetParameter(1, maxPosition);
        primaryFitFunc->SetParameter(2, 0.03);
        primaryFitFunc->SetParameter(3, 1.00);
        primaryFitFunc->SetParameter(4, 0.00);
        primaryFitFunc->SetParameter(5, 0.70);
        primaryFitFunc->SetParameter(6, 1.60);

        // Ограничения для первой аппроксимации
        primaryFitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        primaryFitFunc->SetParLimits(1, 0.85, 1.05);
        primaryFitFunc->SetParLimits(2, 0.0085, 0.085);
        primaryFitFunc->SetParLimits(3, 0.33, 3.00);
        primaryFitFunc->SetParLimits(5, 0.60, 0.95);
        primaryFitFunc->SetParLimits(6, 0.60, 4.00);

        hist->Fit(primaryFitFunc, "RQ");

        // Копируем параметры
        for (int f = 0; f < 7; f++) {
            fitFunc->SetParameter(f, primaryFitFunc->GetParameter(f));
        }

        if (getW(j) > W_RESOLUTION_THRESHOLD) {
            // --- Подгонка с Δ-пиком (второй пик виден) ---
            fitFunc->SetParameter(7, 1.00);
            fitFunc->FixParameter(8, deltaPeak);
            fitFunc->SetParameter(9, 0.03);

            fitFunc->SetParLimits(0, 0.7 * primaryFitFunc->GetParameter(0), 1.3 * primaryFitFunc->GetParameter(0));
            fitFunc->SetParLimits(1, 0.9 * primaryFitFunc->GetParameter(1), 1.1 * primaryFitFunc->GetParameter(1));
            fitFunc->SetParLimits(2, 0.9 * primaryFitFunc->GetParameter(2), 1.1 * primaryFitFunc->GetParameter(2));
            fitFunc->SetParLimits(3, 0.9 * primaryFitFunc->GetParameter(3), 1.1 * primaryFitFunc->GetParameter(3));
            fitFunc->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(deltaPeak)));
            fitFunc->SetParLimits(9, 0.01, 0.1);

            hist->Fit(fitFunc, "RQ");
        } else {
            // --- Δ-пик не виден: задаём значения по умолчанию ---
            fitFunc->FixParameter(7, 0.0);   // amp
            fitFunc->FixParameter(8, 0.0);  // mean
            fitFunc->FixParameter(9, 0.0);  // sigma
        }

        delete primaryFitFunc;
        return fitFunc;
    }


    double getUpEdge(int i, int j, int k, int l) const {
        if (getW(j) > W_RESOLUTION_THRESHOLD) return interpolate(MM_UP_EDGES, j);
        return getUpShortEdge(i, j, k, l);
    }

    double getUpShortEdge(int i, int j, int k, int l) const { return interpolate(MM_UP_SHORT_EDGES, j); }

    double getDownEdge(int i, int j, int k, int l) const { 
        if (getW(j) > W_RESOLUTION_THRESHOLD) return 0.80;
        return getDownShortEdge(i, j, k, l);
    }

    double getDownShortEdge(int i, int j, int k, int l) const { return 0.80; }

    double getDeltaPeak(int i, int j, int k, int l) const { return interpolate(DELTA_PEAK_POSITIONS, j); }

    void getNeutronPeak(TH1F* hist, double& maxPosition, double& maxValue) {
        int bin_min = hist->FindBin(0.86);
        int bin_max = hist->FindBin(1.00);

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

    void getMaxBin(TH1F* hist, double& maxPosition, double& maxValue) const {
        int maxBin = hist->GetMaximumBin();
        maxValue = hist->GetBinContent(maxBin);
        maxPosition = hist->GetBinCenter(maxBin);
    }

    double getNeutronIntegral(TH1F* hist) {
        int bin_min = hist->FindBin(0.86);
        int bin_max = hist->FindBin(1.00);

        double integral = 0.0;

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            integral += hist->GetBinContent(bin);
        }
        integral /= hist->GetBinWidth(1);

        return integral;
    }

private:
    static constexpr double W_RESOLUTION_THRESHOLD = 1.35;

    static constexpr double MM_UP_EDGES[36] = {
        0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.14, 1.18, 1.23, 1.24, 1.24, 1.28, 1.30,
        1.32, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34,
        1.34, 1.34, 1.34, 1.34
    };

    static constexpr double MM_UP_SHORT_EDGES[36] = {
        0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.12, 1.12, 1.10, 1.10, 1.10, 1.10, 1.10,
        1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10,
        1.10, 1.10, 1.10, 1.10
    };

    static constexpr double DELTA_PEAK_POSITIONS[36] = {
        1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.088, 1.088, 1.134, 1.134, 1.169, 1.169, 1.201, 1.201, 
        1.223, 1.223, 1.239, 1.239, 1.230, 1.230, 1.227, 1.227, 1.233, 1.233, 1.216, 1.216, 1.220, 1.220, 1.220, 1.220, 1.220, 1.220
    };

    double interpolate(const double* array, int j) const {
        const int n = 36;
        double x = static_cast<double>(j) * (n - 1) / (NUMBER_W - 1);
        int i = static_cast<int>(x);
        double frac = x - i;

        if (i >= n - 1) return array[n - 1];
        return array[i] + frac * (array[i + 1] - array[i]);
    }

    double getW(int j) const {
        return W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W;
    }
};
