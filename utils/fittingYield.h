#pragma once

#include "TGraphErrors.h"
#include "TF1.h"
#include <cmath>

class FittingYield {
public:
    TF1* fit(TGraphErrors* graph, int i, int j, int k) {
        std::string strFitFunc = "[0] + [1] * cos(x) + [2] * cos(2 * x)";
        TF1* fitFunc = new TF1("Wave Fit Function", strFitFunc.c_str(), 0.0, 2 * M_PI);

        double c, frtAmp, sndAmp;
        getInitialParams(graph, i, j, k, c, frtAmp, sndAmp);

        fitFunc->SetParameter(0, c);
        fitFunc->SetParameter(1, frtAmp);
        fitFunc->SetParameter(2, sndAmp);

        graph->Fit(fitFunc, "RQ");
        return fitFunc;
    }

private:
    void getInitialParams(TGraphErrors* graph, int i, int j, int k,
                          double& c, double& frtAmp, double& sndAmp) {
        if (!graph) {
            c = 0.0;
            frtAmp = 0.0;
            sndAmp = 0.0;
            return;
        }

        int nPoints = graph->GetN();
        if (nPoints == 0) {
            c = 0.0;
            frtAmp = 0.0;
            sndAmp = 0.0;
            return;
        }

        double x, y;
        double maxY = -1e9;
        double minY = 1e9;

        // Поиск минимума и максимума
        for (int ip = 0; ip < nPoints; ++ip) {
            graph->GetPoint(ip, x, y);
            if (y > maxY) maxY = y;
            if (y < minY) minY = y;
        }

        // Задание параметров
        c = 0.5 * (maxY + minY);
        sndAmp = -0.5 * (maxY - minY);
        frtAmp = 0.0;  // начальное приближение для cos(x)
    }
};
