#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

#include "MM_project_utils.h"
#include "rootFitAnalysisStep.h"

class DrawHistStep : public RootFitAnalysisStep {
public:
    DrawHistStep(const std::string& stepName,
                 const std::map<std::string, std::string>& inputFileNames,
                 const std::string& outputDirName = "img/")
        : RootFitAnalysisStep(stepName, inputFileNames),
        outputDirName(outputDirName) {}

    virtual ~DrawHistStep() = default;

protected:
    std::string outputDirName;
    size_t eventsProcessed = 0; ///< Сколько событий уже обработано
    size_t eventsSkipped = 0;
    size_t printEvery = 2000;   ///< Частота логирования прогресса (по умолчанию)

    bool initialize() override {
        if (!RootFitAnalysisStep::initialize()) {
            return false;
        }

        if (!inputFiles.count("hist")) {
            log("Input file with role 'hist' not found.", LogLevel::Error);
            return false;
        }

        return true;
    }

    void processEvent() override {
        std::string name = makeCellName(index_q2, index_w, index_cos_theta, index_phi);

        TH1F* hist = nullptr;
        inputFiles["hist"]->GetObject(name.c_str(), hist);

        if (!hist) {
            ++eventsSkipped;
            return;
        }

        TCanvas* canvas = new TCanvas("canvas", "Histogram", 1200, 900);

        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        TF1* fitFunc = new TF1("FitFunc", strFitFunc.c_str(), downEdge , upEdge);
        TF1* subFunc = new TF1("SubFunc", strSubFunc.c_str(), downEdge, upEdge);
        setFunctionParams(fitFunc, subFunc);

        canvas->Divide(1, 1);
        canvas->cd(1);

        hist->Draw();
        fitFunc->Draw("same");
        subFunc->SetLineColor(kBlue);
        subFunc->Draw("same");

        canvas->Update();

        std::string figName = outputDirName + makeCellName(index_q2, index_w, index_cos_theta, index_phi) + ".png";
        canvas->SaveAs(figName.c_str());

        delete fitFunc;
        delete subFunc;
        delete canvas;
        delete hist;
    }

    bool check() override {
        if (ampGauss < 0.0) return false;
        return true;
    }

    void logProgress() override {
        ++eventsProcessed;
        if (eventsProcessed % printEvery == 0 || eventsProcessed == totalEvents) {
            double percent = 100.0 * (eventsProcessed - eventsSkipped) / totalEvents;
            log("Processed " + std::to_string(eventsProcessed - eventsSkipped) +
                " histograms (" + std::to_string(static_cast<int>(percent)) + "%)", LogLevel::Info);
        }
    }

private:
    void setFunctionParams(TF1* fitFunc, TF1* subFunc) {
        fitFunc->SetParameter(0, ampGauss);
        fitFunc->SetParameter(1, meanGauss);
        fitFunc->SetParameter(2, stdDevGauss);
        fitFunc->SetParameter(3, asymGauss);
        fitFunc->SetParameter(4, paramA);
        fitFunc->SetParameter(5, paramB);
        fitFunc->SetParameter(6, paramC);
        fitFunc->SetParameter(7, ampGaussDelta);
        fitFunc->SetParameter(8, meanGaussDelta);
        fitFunc->SetParameter(9, stdDevGaussDelta);

        subFunc->SetParameter(4, paramA);
        subFunc->SetParameter(5, paramB);
        subFunc->SetParameter(6, paramC);
    }
};