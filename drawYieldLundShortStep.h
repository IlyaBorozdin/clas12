#pragma once

#include "TTree.h"
#include "TGraphErrors.h"

#include "source/ManageClasses/rootListObjectAnalysisStep.h"
#include "MM_project_utils.h"

class DrawYieldLundShortStep : public RootListObjectAnalysisStep<TTree> {
public:
    DrawYieldLundShortStep(const std::string& stepName,
                           const std::string& inputFileName)
        : RootListObjectAnalysisStep<TTree>(stepName, { {"main", inputFileName} }, cellGenerator) {}

    ~DrawYieldLundShortStep() override = default;

    std::vector<std::vector<std::vector<TGraphErrors*>>> getGraphs() const {
        return graphs;
    }

protected:
    size_t processedTrees = 0;
    size_t skippedTrees = 0;
    size_t printEvery = 2000;

    std::vector<std::vector<std::vector<TGraphErrors*>>> graphs;

    bool initialize() override {
        if (!RootListObjectAnalysisStep<TTree>::initialize()) return false;

        // Создаём графики
        graphs.resize(NUMBER_Q2);
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            graphs[iq2].resize(NUMBER_W);
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                graphs[iq2][iw].resize(NUMBER_COS_THETA, nullptr);
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    graphs[iq2][iw][ict] = new TGraphErrors(NUMBER_PHI);
                    graphs[iq2][iw][ict]->SetTitle(generateTitleForYieldPlot(iq2, iw, ict).c_str());
                    graphs[iq2][iw][ict]->SetLineWidth(2);
                    graphs[iq2][iw][ict]->SetLineStyle(1);
                    graphs[iq2][iw][ict]->SetLineColor(kGreen);

                    for (int iphi = 0; iphi < NUMBER_PHI; ++iphi) {
                        double phiCenter = (iphi + 0.5) * (2 * M_PI / NUMBER_PHI);
                        double phiError = (M_PI / NUMBER_PHI);
                        graphs[iq2][iw][ict]->SetPoint(iphi, phiCenter, 0.0);  // Y = 0
                    }
                }
            }
        }

        return true;
    }

    bool finalize() override {
        log("Processed trees: " + std::to_string(processedTrees), LogLevel::Info);
        log("Skipped trees due to missing branch or bad name: " + std::to_string(skippedTrees), LogLevel::Info);

        return RootListObjectAnalysisStep<TTree>::finalize();
    }

    void processObject(TTree* tree) override {
        // Извлекаем координаты ячейки
        int index_q2, index_w, index_cos_theta, index_phi;
        if (!parseCellName(tree->GetName(), index_q2, index_w, index_cos_theta, index_phi)) {
            ++skippedTrees;
            return;
        }

        Long64_t nEntries = tree->GetEntries();

        // Вычисление выхода реакции и ошибки
        double yPoint = nEntries;
        double yErr = std::sqrt(nEntries);

        double phiCenter = (index_phi + 0.5) * (2 * M_PI / NUMBER_PHI); // центр бина по φ
        double phiError = (M_PI / NUMBER_PHI); // половина ширины бина

        // Добавление точки в соответствующий график
        graphs[index_q2][index_w][index_cos_theta]->SetPoint(index_phi, phiCenter, yPoint);
        graphs[index_q2][index_w][index_cos_theta]->SetPointError(index_phi, phiError, yErr);

    }

    void logProgress() override {
        ++processedTrees;
        if (processedTrees % printEvery == 0 || processedTrees == totalObjects) {
            log("Processing tree " + std::to_string(processedTrees) + " / " + std::to_string(totalObjects), LogLevel::Info);
        }
    }
};