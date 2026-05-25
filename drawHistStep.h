#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"

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

    size_t eventsRejectedCondition4 = 0;
    size_t eventsRejectedCondition5 = 0;
    size_t eventsRejectedCondition6 = 0;
    size_t eventsRejectedCondition7 = 0;


    bool initialize() override {
        if (!RootFitAnalysisStep::initialize()) {
            return false;
        }

        if (!inputFiles.count("hist")) {
            log("Input file with role 'hist' not found.", LogLevel::Error);
            return false;
        }

        // --- Глобальные настройки стиля для отрисовки PNG ---
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "t");
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleSize(0.04, "t");
        gStyle->SetTitleSize(0.04, "XYZ");
        gStyle->SetLabelSize(0.035, "XYZ");
        // Чтобы русские буквы в заголовках осей не обрезались
        gStyle->SetTitleOffset(1.2, "Y");

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

        // Формулы фитирования
        std::string strSubFunc = "TMath::Max(0., [4]*(x - [5])*(x - [6]))";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        TF1* fitFunc = new TF1("FitFunc", strFitFunc.c_str(), downEdge , upEdge);
        TF1* subFunc = new TF1("SubFunc", strSubFunc.c_str(), downEdge, upEdge);
        setFunctionParams(fitFunc, subFunc);

        canvas->SetGrid();
        canvas->cd();

        // Рисуем гистограмму. Т.к. первая панель уже в списке функций hist, она отрисуется сама
        hist->Draw(); 
        
        fitFunc->Draw("same");
        subFunc->SetLineColor(kBlue);
        subFunc->SetLineStyle(2); // Пунктир для фона выглядит нагляднее
        subFunc->Draw("same");

        // --- СОЗДАНИЕ ВТОРОЙ ПАНЕЛИ (Yield) ---
        // Координаты первой панели были (0.65, 0.65, 0.90, 0.90). 
        // Ставим новую панель чуть ниже: от 0.52 до 0.63 по вертикали.
        TPaveText* yieldInfo = new TPaveText(0.65, 0.52, 0.90, 0.63, "NDC");
        yieldInfo->SetBorderSize(1);
        yieldInfo->SetFillColor(0);
        yieldInfo->SetTextAlign(12); // Слева по центру
        yieldInfo->SetTextFont(42);

        yieldInfo->SetTextSize(0.035); 
        yieldInfo->SetMargin(0.05); // Небольшой внутренний отступ для красоты

        // Рассчитываем относительную погрешность
        double relError = (yield != 0) ? (eYield / yield) : 0.0;

        // --- В блоке создания Yield панели ---
        if (yield > 0) {
            double relError = eYield / yield;
            yieldInfo->AddText(Form("Y = %.2f #pm %.2f", yield, eYield));
            yieldInfo->AddText(Form("#deltaY / Y = %.3f", relError));
        } else {
            // Для пустых бинов показываем только верхний предел
            yieldInfo->AddText("Y = 0.00");
            yieldInfo->AddText(Form("up limit: %.2f", eYield)); 
        }
        
        yieldInfo->Draw("same"); // Отрисовываем поверх

        canvas->Update();

        std::string figName = outputDirName + name + ".png";
        canvas->SaveAs(figName.c_str());

        // Очистка
        delete yieldInfo;
        delete fitFunc;
        delete subFunc;
        delete canvas;
    }

    bool check() override {
        bool passed = true;

        /*
        double valueAtPeak = std::max(0.0, paramA * (meanGauss - paramB) * (meanGauss - paramC)) + ampGauss;
        if (valueAtPeak > 1.25 * neutronValue) {
            ++eventsRejectedCondition4;
            passed = false;
        }

        if (valueAtPeak < 0.75 * neutronValue) {
            ++eventsRejectedCondition5;
            passed = false;
        }
        */

        /*
        if (ampGaussDelta > 1.25 * maxValue) {
            ++eventsRejectedCondition6;
            passed = false;
        }

        if (ampGaussDelta < 0.0) {
            ++eventsRejectedCondition7;
            passed = false;
        }
        */

        return passed;
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
        log("Histograms drawing complete.", LogLevel::Info);
        log("Processed histogramms: " + std::to_string(eventsProcessed), LogLevel::Info);
        log("Skipped histogramms due to empty or bad name: " + std::to_string(eventsSkipped), LogLevel::Info);
        /*
        log("Rejected by Condition 4 (peak value > 1.25 * neutron): " + std::to_string(eventsRejectedCondition4), LogLevel::Info);
        log("Rejected by Condition 5 (peak value < 0.75 * neutron): " + std::to_string(eventsRejectedCondition5), LogLevel::Info);
        log("Rejected by Condition 6 (delta peak too high): " + std::to_string(eventsRejectedCondition6), LogLevel::Info);
        log("Rejected by Condition 7 (delta peak negative): " + std::to_string(eventsRejectedCondition7), LogLevel::Info);
        */

        return RootFitAnalysisStep::finalize();
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