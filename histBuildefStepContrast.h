#pragma once

#include <cmath>

#include "TH1F.h"
#include "TTree.h"

#include "source/ManageClasses/rootListObjectAnalysisStep.h"
#include "MM_project_const.h"
#include "MM_project_utils.h"
#include "utils/binOptimizer.h"

/**
 * @brief Шаг анализа, строящий гистограммы по недостающей массе (MM) 
 *        для каждого дерева в 4D-сетке (Q², W, cosθ, φ) с использованием оптимизации контрастности.
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

    bool finalize() override {
        log("Histograms building complete.", LogLevel::Info);
        log("Processed trees: " + std::to_string(processedTrees), LogLevel::Info);
        log("Skipped trees due to missing branch or bad name: " + std::to_string(skippedTrees), LogLevel::Info);

        return RootListObjectAnalysisStep<TTree>::finalize();
    }

    void processObject(TTree* tree) override {
        // Извлекаем координаты ячейки
        int i, j, k, l;
        if (!parseCellName(tree->GetName(), i, j, k, l)) {
            ++skippedTrees;
            return;
        }

        // Читаем данные из дерева
        double MM;
        tree->SetBranchAddress("MM", &MM);
        std::vector<double> values;
        Long64_t nEntries = tree->GetEntries();
        values.reserve(nEntries);

        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            values.push_back(MM);
        }

        // Оптимизация числа бинов
        BinOptimizer optimizer(values, 5.0, 0.5);
        int optimalBins = optimizer.findOptimalBinCount();

        std::string name = makeCellName(i, j, k, l);
        std::string title = generateTitle(i, j, k, l);

        // Построение гистограммы
        double xMin = optimizer.getXMin();
        double xMax = optimizer.getXMax();
        TH1F* hist = new TH1F(name.c_str(), title.c_str(), optimalBins, xMin, xMax);

        for (double val : values) {
            if (val >= xMin && val <= xMax)
                hist->Fill(val);
        }

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
    /**
     * @brief Генерация подписи гистограммы по координатам ячейки
     */
    std::string generateTitle(int i, int j, int k, int l) const {
        double q2Min = STEPS_Q2[2 * i];
        double q2Max = STEPS_Q2[2 * i + 1];

        double wMin = W_MIN + STEP_W * j;
        double wMax = W_MIN + STEP_W * (j + 1);

        double cosThetaMin = -1.0 + STEP_COS_THETA * k;
        double cosThetaMax = -1.0 + STEP_COS_THETA * (k + 1);

        double phiMinDeg = STEP_PHI * l * 180.0 / M_PI;
        double phiMaxDeg = STEP_PHI * (l + 1) * 180.0 / M_PI;

        return "MM: Q^{2} = " + formatFloat(q2Min, 1) + "-" + formatFloat(q2Max, 1) +
               " GeV^{2}, W = " + formatFloat(wMin, 2) + "-" + formatFloat(wMax, 2) + " GeV, " +
               "cos(Theta) = (" + formatFloat(cosThetaMin, 1) + ")-(" + formatFloat(cosThetaMax, 1) + "), " +
               "Phi = " + formatFloat(phiMinDeg, 0) + "^{o}-" + formatFloat(phiMaxDeg, 0) + "^{o}; MM, GeV";
    }
};
