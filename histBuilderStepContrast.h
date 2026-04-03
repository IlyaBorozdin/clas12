#pragma once

#include <cmath>
#include <vector>

#include "TH1F.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TList.h"
#include "TStyle.h"

#include "source/ManageClasses/rootListObjectAnalysisStep.h"
#include "MM_project_const.h"
#include "MM_project_utils.h"
#include "utils/binOptimizer.h"

/**
 * @brief Analysis step for building Missing Mass (MM) histograms 
 * with automated Y-axis scaling to accommodate info panels.
 */
class HistBuilderStepContrast : public RootListObjectAnalysisStep<TTree> {
public:
    HistBuilderStepContrast(const std::string& stepName,
                            const std::string& inputFileName,
                            const std::string& outputFileName)
        : RootListObjectAnalysisStep<TTree>(stepName, { {"main", inputFileName} }, cellGenerator, outputFileName) {}

    ~HistBuilderStepContrast() override = default;

protected:
    size_t processedTrees = 0;
    size_t skippedTrees = 0;
    size_t printEvery = 2000;

    bool initialize() override {
        if (!RootListObjectAnalysisStep<TTree>::initialize()) {
            return false;
        }

        // Global Style Settings
        gStyle->SetTitleFont(42, "t");
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTextFont(42);
        
        gStyle->SetTitleSize(0.045, "t");
        gStyle->SetTitleSize(0.045, "XYZ");
        gStyle->SetLabelSize(0.035, "XYZ");

        return true;
    }

    bool finalize() override {
        log("Histograms building complete.", LogLevel::Info);
        log("Processed trees: " + std::to_string(processedTrees), LogLevel::Info);
        log("Skipped trees: " + std::to_string(skippedTrees), LogLevel::Info);

        return RootListObjectAnalysisStep<TTree>::finalize();
    }

    void processObject(TTree* tree) override {
        int i, j, k, l;
        if (!parseCellName(tree->GetName(), i, j, k, l)) {
            ++skippedTrees;
            return;
        }

        double MM = 0.0;
        double weight = 1.0;
        tree->SetBranchAddress("MM", &MM);
        TBranch* weightBranch = tree->GetBranch("weight");
        if (weightBranch) tree->SetBranchAddress("weight", &weight);

        std::vector<WeightedValue> data;
        Long64_t nEntries = tree->GetEntries();
        data.reserve(nEntries);
        
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            data.emplace_back(MM, weightBranch ? weight : 1.0);
        }

        BinOptimizer optimizer(data, 5.0, 0.5);
        int optimalBins = optimizer.findOptimalBinCount();
        double xMin = optimizer.getXMin();
        double xMax = optimizer.getXMax();

        std::string name = makeCellName(i, j, k, l);

        // --- 1. Create Histogram with English titles ---
        std::string fullTitle = ";MM(e^{-}p#rightarrow e^{-}#pi^{+}X), GeV;Counts";
        TH1F* hist = new TH1F(name.c_str(), fullTitle.c_str(), optimalBins, xMin, xMax);
        hist->SetStats(0); 

        for (const auto& d : data) {
            if (d.value >= xMin && d.value <= xMax)
                hist->Fill(d.value, d.weight);
        }

        // --- 2. Dynamic Y-axis Scaling ---
        // Добавляем 35% запаса сверху, чтобы TPaveText ничего не перекрывал
        double maxContent = hist->GetMaximum();
        hist->SetMaximum(maxContent * 1.35);
        hist->SetMinimum(0);

        // --- 3. Info Pave Panel ---
        // Координаты NDC (от 0 до 1), панель в верхнем правом углу
        TPaveText* info = new TPaveText(0.62, 0.65, 0.88, 0.88, "NDC");
        info->SetBorderSize(1);
        info->SetFillColor(0);
        info->SetTextAlign(12); // Left-aligned, vertically centered
        info->SetTextFont(42);
        info->SetTextSize(0.035); // Фиксированный размер текста

        addKinematicInfo(info, i, j, k, l, static_cast<Long64_t>(hist->GetEntries()));

        // Добавляем в список функций, чтобы рисовалось автоматически при Draw()
        hist->GetListOfFunctions()->Add(info);

        outputFile->cd();
        hist->Write();
    }

    void logProgress() override {
        ++processedTrees;
        if (processedTrees % printEvery == 0 || processedTrees == totalObjects) {
            log("Processing tree " + std::to_string(processedTrees) + " / " + std::to_string(totalObjects), LogLevel::Info);
        }
    }

private:
    void addKinematicInfo(TPaveText* pave, int i, int j, int k, int l, Long64_t entries) const {
        double q2Min = STEPS_Q2[2 * i], q2Max = STEPS_Q2[2 * i + 1];
        double wMin = W_MIN + STEP_W * j, wMax = W_MIN + STEP_W * (j + 1);
        double cosMin = -1.0 + STEP_COS_THETA * k, cosMax = -1.0 + STEP_COS_THETA * (k + 1);
        double phiMin = STEP_PHI * l * 180.0 / M_PI, phiMax = STEP_PHI * (l + 1) * 180.0 / M_PI;

        pave->AddText(Form("Entries: %lld", entries)); 
        pave->AddText(Form("Q^{2}: %.2f-%.2f GeV^{2}", q2Min, q2Max));
        pave->AddText(Form("W: %.3f-%.3f GeV", wMin, wMax));
        pave->AddText(Form("cos(#theta^{*}): %.2f-%.2f", cosMin, cosMax));
        pave->AddText(Form("#varphi^{*}: %.0f^{o}-%.0f^{o}", phiMin, phiMax));
    }
};