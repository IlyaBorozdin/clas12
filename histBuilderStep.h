#pragma once

#include <cmath>

#include "TH1F.h"
#include "TTree.h"

// #include "source/ManageClasses/rootListDataAnalysisStep.h"
#include "source/ManageClasses/rootListObjectAnalysisStep.h"
#include "MM_project_const.h"
#include "MM_project_utils.h"

/**
 * @brief Шаг анализа, строящий гистограммы по недостающей массе (MM) 
 *        для каждого дерева в 4D-сетке (Q², W, cosθ, φ).
 *
 * Класс предполагает, что входной ROOT-файл содержит множество деревьев,
 * имена которых кодируют позицию в сетке (i,j,k,l). Каждое дерево содержит
 * ветку "MM", и по ней строится одномерная гистограмма, записываемая в выходной файл.
 */
class HistBuilderStep : public RootListObjectAnalysisStep<TTree> {
public:
    HistBuilderStep(const std::string& stepName,
                    const std::string& inputFileName,
                    const std::string& outputFileName)
        : RootListObjectAnalysisStep<TTree>(stepName, { {"main", inputFileName} }, cellGenerator, outputFileName) {}

    virtual ~HistBuilderStep() = default;

protected:
    size_t processedTrees = 0;   ///< Счётчик обработанных деревьев
    size_t skippedTrees = 0;     ///< Счётчик пропущенных деревьев
    size_t printEvery = 2000;    ///< Частота логирования прогресса

    /**
     * @brief Завершение шага: закрытие файла и отчёт
     */
    bool finalize() override {
        log("Histograms building complete.", LogLevel::Info);
        log("Processed trees: " + std::to_string(processedTrees), LogLevel::Info);
        log("Skipped trees due to missing branch or bad name: " + std::to_string(skippedTrees), LogLevel::Info);

        return RootListObjectAnalysisStep<TTree>::finalize();
    }

    /**
     * @brief Обработка одного дерева: строим и сохраняем гистограмму MM
     */
    void processObject(TTree* tree) override {
        // Извлекаем (i, j, k, l) из имени
        int i, j, k, l;
        if (!parseCellName(tree->GetName(), i, j, k, l)) {
            ++skippedTrees;
            return;
        }

        Long64_t nEntries = tree->GetEntries();
        int optimalBins = calculateOptimalBins(nEntries);
        std::string name = makeCellName(i, j, k, l);
        std::string title = generateTitle(i, j, k, l);

        TH1F* hist = new TH1F(name.c_str(), title.c_str(), optimalBins, MM_MIN, MM_MAX);

        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            hist->Fill(MM);
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
    double MM;

    /**
     * @brief Выбор оптимального количества бинов по числу событий
     */
    int calculateOptimalBins(size_t N) {
        const int B_MIN = 20;
        const int B_MAX = 80;
        const double N0 = 15000.0;

        double growthFactor = 1.0 - std::exp(-static_cast<double>(N) / N0);
        int bins = static_cast<int>(B_MIN + (B_MAX - B_MIN) * growthFactor);
        return std::max(B_MIN, bins);
    }

    /**
     * @brief Генерация подписи гистограммы по координатам ячейки
     */
    std::string generateTitle(int i, int j, int k, int l) const {
        // Границы по Q²
        double q2Min = STEPS_Q2[2 * i];
        double q2Max = STEPS_Q2[2 * i + 1];

        // Границы по W
        double wMin = W_MIN + STEP_W * j;
        double wMax = W_MIN + STEP_W * (j + 1);

        // Границы по cos(θ)
        double cosThetaMin = -1.0 + STEP_COS_THETA * k;
        double cosThetaMax = -1.0 + STEP_COS_THETA * (k + 1);

        // Границы по φ (в градусах)
        double phiMinDeg = STEP_PHI * l * 180.0 / M_PI;
        double phiMaxDeg = STEP_PHI * (l + 1) * 180.0 / M_PI;

        return "MM: Q^{2} = " + formatFloat(q2Min, 1) + "-" + formatFloat(q2Max, 1) +
               " GeV^{2}, W = " + formatFloat(wMin, 2) + "-" + formatFloat(wMax, 2) + " GeV, " +
               "cos(Theta) = (" + formatFloat(cosThetaMin, 1) + ")-(" + formatFloat(cosThetaMax, 1) + "), " +
               "Phi = " + formatFloat(phiMinDeg, 0) + "^{o}-" + formatFloat(phiMaxDeg, 0) + "^{o}; MM, GeV";
    }
};
