#pragma once

#include "TGraphErrors.h"
#include "TF1.h"

class FittingYield {
public:
    TF1* fit(TGraphErrors* graph, int i, int j, int k) {
        std::string strFitFunc = "[0] + [1] * cos(x) + [2] * cos(2 * x)";
        TF1* fitFunc = new TF1("Wave Fit Function", strFitFunc.c_str(), 0.0, 2 * M_PI);

        double sndAmp;
        getSecondAmp(graph, i, j, k, sndAmp);

        fitFunc->SetParameter(0, 0.0);
        fitFunc->SetParameter(1, 0.0);
        fitFunc->SetParameter(2, sndAmp);

        graph->Fit(fitFunc, "RQ");
        return fitFunc;
    }

    void getSecondAmp(TGraphErrors* graph, int i, int j, int k, double& sndAmp) {
        if (!graph) {
            sndAmp = 0.0;
            return;
        }

        int nPoints = graph->GetN();
        if (nPoints == 0) {
            sndAmp = 0.0;
            return;
        }

        double maxY = 0.0;
        double x, y;

        // Поиск максимального значения Y
        for (int ip = 0; ip < nPoints; ++ip) {
            graph->GetPoint(ip, x, y);
            if (y > maxY) {
                maxY = y;
            }
        }

        sndAmp = -maxY;
        }
};