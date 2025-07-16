#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2I.h"
#include "TFile.h"

#include "MM_project_utils.h"
#include "MM_project_const.h"
#include "rootFitAnalysisStep.h"

class DrawHistStepDebug : public RootFitAnalysisStep {
public:
    DrawHistStepDebug(const std::string& stepName,
                      const std::string& inputFileName,
                      const std::string& debugFileName)
        : RootFitAnalysisStep(stepName,{ { "main", inputFileName } }, debugFileName) {}

    virtual ~DrawHistStepDebug() = default;

protected:
    size_t eventsProcessed = 0; ///< Сколько событий уже обработано
    size_t eventsSkipped = 0;
    size_t printEvery = 2000;   ///< Частота логирования прогресса (по умолчанию)

    TH2I* condition4Maps[NUMBER_Q2][NUMBER_W] = {};

    bool initialize() override {
        if (!RootFitAnalysisStep::initialize()) return false;

        // Создаём все 72 карты
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                std::string name = "condsMap_q2_" + std::to_string(iq2 + 1) + "_w_" + std::to_string(iw + 1);
                condition4Maps[iq2][iw] = new TH2I(name.c_str(), name.c_str(),
                                                   NUMBER_PHI, 0, 2 * M_PI,
                                                   NUMBER_COS_THETA, -1.0, 1.0);
                condition4Maps[iq2][iw]->Reset();
            }
        }

        return true;
    }

    void processEvent() override {
        // Проверка всех условий (аналог check())
        bool passed = true;

        // Защита от отрицательного значения в параболе
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

        // Проверка на корректные индексы
        if (index_q2 >= NUMBER_Q2 || index_w >= NUMBER_W ||
            index_phi >= NUMBER_PHI || index_cos_theta >= NUMBER_COS_THETA) {
            ++eventsSkipped;
            return;
        }

        // Запись результата в карту
        condition4Maps[index_q2][index_w]->SetBinContent(index_phi + 1, index_cos_theta + 1, passed ? 1 : 0);
    }

    void logProgress() override {
        ++eventsProcessed;
        if (eventsProcessed % printEvery == 0 || eventsProcessed == totalEvents) {
            double percent = 100.0 * (eventsProcessed) / totalEvents;
            log("Processed " + std::to_string(eventsProcessed) +
                " histograms (" + std::to_string(static_cast<int>(percent)) + "%)", LogLevel::Info);
        }
    }

    bool finalize() override {
        log("Debug histogram drawing complete.", LogLevel::Info);
        log("Processed histogramms: " + std::to_string(eventsProcessed), LogLevel::Info);
        log("Skipped histogramms due to empty or bad name: " + std::to_string(eventsSkipped), LogLevel::Info);

        outputFile->cd();
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                if (condition4Maps[iq2][iw]) {
                    condition4Maps[iq2][iw]->Write();
                }
            }
        }

        return RootFitAnalysisStep::finalize();
    }
};
