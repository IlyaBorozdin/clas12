#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "MM_project_utils.h"
#include "genericFitAnalysisStep.h"

class NewDrawHistStep : public GenericFitAnalysisStep {
public:
    NewDrawHistStep(const std::string& stepName,
                    const std::map<std::string, std::string>& inputFileNames,
                    const std::string& outputDirName,
                    std::unique_ptr<FitModel> fitModel)
        : GenericFitAnalysisStep(stepName, inputFileNames, "", std::move(fitModel)),
          outputDirName(outputDirName) {}

protected:
    std::string outputDirName;
    size_t eventsSkipped = 0;

    bool initialize() override {
        if (!GenericFitAnalysisStep::initialize()) return false;

        if (!inputFiles.count("hist")) {
            log("Input file with role 'hist' not found.", LogLevel::Error);
            return false;
        }

        // --- Настройка стиля (как в старом классе) ---
        gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42, "t");
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTitleSize(0.04, "t");
        gStyle->SetTitleSize(0.04, "XYZ");
        gStyle->SetLabelSize(0.035, "XYZ");
        gStyle->SetTitleOffset(1.2, "Y");

        return true;
    }

    void processEvent() override {
        // Формируем имя гистограммы по индексам из дерева
        std::string name = makeCellName(index_q2, index_w, index_cos_theta, index_phi);

        TH1F* hist = nullptr;
        inputFiles["hist"]->GetObject(name.c_str(), hist);

        if (!hist) {
            ++eventsSkipped;
            return;
        }

        TCanvas* canvas = new TCanvas("c", "Fit Display", 1200, 900);
        canvas->SetGrid();
        canvas->cd();

        // 1. Создаем функции через полиморфную модель
        TF1* fitFunc = model->createTF1("FitFunc", downEdge, upEdge);
        TF1* subFunc = model->getBackgroundFormula("SubFunc", downEdge, upEdge);

        // 2. Устанавливаем параметры (модель берет их из своих внутренних буферов, 
        // которые заполнились при GetEntry)
        model->setParams(fitFunc, subFunc);

        // 3. Отрисовка
        hist->SetTitle(name.c_str());
        hist->Draw();
        
        fitFunc->SetLineColor(kRed);
        fitFunc->SetLineWidth(3);
        fitFunc->Draw("same");

        if (subFunc) {
            subFunc->SetLineColor(kBlue);
            subFunc->SetLineStyle(2);
            subFunc->Draw("same");
        }

        // 4. Инфо-панель (универсальная часть)
        TPaveText* yieldInfo = drawInfoPave();
        yieldInfo->Draw("same");

        canvas->Update();
        std::string figPath = outputDirName + name + ".png";
        canvas->SaveAs(figPath.c_str());

        // Чистка
        delete yieldInfo;
        delete fitFunc;
        delete subFunc;
        delete canvas;
    }

    bool finalize() override {
        log("Histograms drawing complete.", LogLevel::Info);
        log("Skipped histograms: " + std::to_string(eventsSkipped), LogLevel::Info);
        return GenericFitAnalysisStep::finalize();
    }

private:
    TPaveText* drawInfoPave() {
        double y2 = (yield > 0) ? 0.68 : 0.63; // 3 строки vs 2 строки
        TPaveText* pt = new TPaveText(0.65, 0.52, 0.90, y2, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.035);
        pt->SetMargin(0.05); 

        // Используем yield и eYield из GenericFitAnalysisStep
        if (yield > 0) {
            double relError = eYield / yield;
            pt->AddText(Form("Y = %.2f #pm %.2f", yield, eYield));
            pt->AddText(Form("#deltaY / Y = %.3f", relError));
            TText *t = pt->AddText(Form("#chi^{2} / NDF = %.2f", chi2ndf));
            if (chi2ndf > 3.0) t->SetTextColor(kRed);
        } else {
            pt->AddText("Y = 0.00");
            pt->AddText(Form("Up limit: %.2f", eYield));
        }

        return pt;
    }
};