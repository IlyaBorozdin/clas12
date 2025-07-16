#pragma once

#include "TGraphErrors.h"
#include "TCanvas.h"

#include "MM_project_utils.h"
#include "MM_project_const.h"
#include "rootFitAnalysisStep.h"

class DrawYieldStep : public RootFitAnalysisStep {
public:
    DrawYieldStep(const std::string& stepName,
                      const std::string& inputFileName,
                      const std::string& outputDirName = "img/")
        : RootFitAnalysisStep(stepName,{ { "main", inputFileName } }),
        outputDirName(outputDirName) {}

    virtual ~DrawYieldStep() = default;

protected:
    std::string outputDirName;
    size_t eventsProcessed = 0;
    size_t printEvery = 2000;

    std::vector<std::vector<std::vector<TGraphErrors*>>> graphs;
    std::vector<std::vector<TCanvas*>> canvases;


    bool initialize() override {
        if (!RootFitAnalysisStep::initialize()) return false;

        // Создаём графики и канвасы
        graphs.resize(NUMBER_Q2);
        canvases.resize(NUMBER_Q2);
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            graphs[iq2].resize(NUMBER_W);
            canvases[iq2].resize(NUMBER_W, nullptr);
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                graphs[iq2][iw].resize(NUMBER_COS_THETA, nullptr);
                std::string cname = "canvas_q2_" + std::to_string(iq2 + 1) + "_w_" + std::to_string(iw + 1);
                canvases[iq2][iw] = new TCanvas(cname.c_str(), cname.c_str(), 1600, 1200);
                canvases[iq2][iw]->Divide(4, 3);

                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    graphs[iq2][iw][ict] = new TGraphErrors(NUMBER_PHI);
                    graphs[iq2][iw][ict]->SetTitle(generateTitleForYieldPlot(iq2, iw, ict).c_str());
                    graphs[iq2][iw][ict]->SetLineWidth(2);
                    graphs[iq2][iw][ict]->SetLineStyle(1);
                    graphs[iq2][iw][ict]->SetLineColor(kBlue);

                    for (int iphi = 0; iphi < NUMBER_PHI; ++iphi) {
                        double phiCenter = (iphi + 0.5) * (2 * M_PI / NUMBER_PHI);
                        double phiError = (M_PI / NUMBER_PHI);
                        graphs[iq2][iw][ict]->SetPoint(iphi, phiCenter, 0.0);  // Y = 0
                        // graphs[iq2][iw][ict]->SetPointError(iphi, phiError, 0.0);
                    }
                }
            }
        }

        return true;
    }

    void processEvent() override {
        // Вычисление выхода реакции и ошибки
        double yPoint = std::sqrt(M_PI / 2) * ampGauss * std::abs(stdDevGauss) * (1 + std::abs(asymGauss)) / width_mm;
        double yErr = yPoint * std::sqrt(
            std::pow(eAmpGauss / ampGauss, 2) +
            std::pow(eStdDevGauss / stdDevGauss, 2) +
            std::pow(eAsymGauss / (1 + asymGauss), 2)
        );

        double phiCenter = (index_phi + 0.5) * (2 * M_PI / NUMBER_PHI); // центр бина по φ
        double phiError = (M_PI / NUMBER_PHI); // половина ширины бина

        // Добавление точки в соответствующий график
        graphs[index_q2][index_w][index_cos_theta]->SetPoint(index_phi, phiCenter, yPoint);
        graphs[index_q2][index_w][index_cos_theta]->SetPointError(index_phi, phiError, yErr);
    }

    void logProgress() override {
        ++eventsProcessed;
        if (eventsProcessed % printEvery == 0 || eventsProcessed == totalEvents) {
            double percent = 100.0 * (eventsProcessed) / totalEvents;
            log("Processed " + std::to_string(eventsProcessed) +
                " histograms (" + std::to_string(static_cast<int>(percent)) + "%)", LogLevel::Info);
        }
    }

    bool check() override {
        bool passed = true;

        double valueAtPeak = std::max(0.0, paramA * (meanGauss - paramB) * (meanGauss - paramC)) + ampGauss;
        if (valueAtPeak > 1.25 * neutronValue) {
            passed = false;
        }

        if (valueAtPeak < 0.75 * neutronValue) {
            passed = false;
        }

        if (ampGaussDelta > 1.25 * maxValue) {
            passed = false;
        }

        if (ampGaussDelta < 0.0) {
            passed = false;
        }

        return passed;
    }

    bool finalize() override {
        log("Yield plotting complete.", LogLevel::Info);
        log("Processed entries: " + std::to_string(eventsProcessed), LogLevel::Info);

        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                // canvas->cd();

                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    canvases[iq2][iw]->cd(ict + 1);
                    graphs[iq2][iw][ict]->Draw("APL");
                    // graphs[iq2][iw][ict]->Write();
                }

                canvases[iq2][iw]->Update();

                std::string figName = "yield_q2_" + std::to_string(iq2 + 1) +
                                    "_w_" + std::to_string(iw + 1) + ".png";
                canvases[iq2][iw]->SaveAs((outputDirName + figName).c_str());
            }
        }

        return RootFitAnalysisStep::finalize();
    }

    std::string generateTitleForYieldPlot(int iq2, int iw, int ict) const {
        // Границы по Q²
        double q2Min = STEPS_Q2[2 * iq2];
        double q2Max = STEPS_Q2[2 * iq2 + 1];

        // Границы по W
        double wMin = W_MIN + STEP_W * iw;
        double wMax = W_MIN + STEP_W * (iw + 1);

        // Границы по cos(θ)
        double cosThetaMin = -1.0 + STEP_COS_THETA * ict;
        double cosThetaMax = -1.0 + STEP_COS_THETA * (ict + 1);

        return "Yield vs #phi: Q^{2} = " + formatFloat(q2Min, 1) + "-" + formatFloat(q2Max, 1) +
            " GeV^{2}, W = " + formatFloat(wMin, 2) + "-" + formatFloat(wMax, 2) + " GeV, " +
            "cos(#theta) = " + formatFloat(cosThetaMin, 1) + " to " + formatFloat(cosThetaMax, 1) + "; #phi, rad; Yield";
    }

};
