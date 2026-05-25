#pragma once

#include "source/manageClasses/rootStep.h"

#include "TChain.h"
#include "TKey.h"
#include "TTree.h"

#include <vector>
#include <string>
#include <iostream>

/**
 * @brief Шаг для объединения всех циклов дерева в один.
 *
 * Из входного файла берёт дерево (по умолчанию "ExpData"),
 * объединяет все его циклы с помощью TChain и сохраняет
 * в новый выходной файл без циклов.
 */
class MergeCyclesStep : public RootStep {
public:
    MergeCyclesStep(const std::string& stepName,
                    const std::string& inputFileName,
                    const std::string& outputFileName)
        : RootStep(stepName, { {"main", inputFileName} }, outputFileName),
          treeName("ExpData") {}

protected:
    std::string treeName;

    bool run() override {
        // Собираем все циклы для указанного дерева
        std::vector<int> cycles;
        TList* keys = inputFile->GetListOfKeys();
        for (TObject* obj : *keys) {
            auto* key = dynamic_cast<TKey*>(obj);
            if (!key) continue;
            if (treeName == std::string(key->GetName())) {
                cycles.push_back(key->GetCycle());
            }
        }

        if (cycles.empty()) {
            log("No cycles for tree '" + treeName + "' found.", LogLevel::Error);
            return false;
        }

        log("Found " + std::to_string(cycles.size()) +
                " cycle(s) for tree '" + treeName + "'.",
            LogLevel::Info);

        // Создаём TChain и добавляем все циклы
        TChain chain(treeName.c_str());
        for (int cyc : cycles) {
            std::string spec = getInputFileName("main") + "/" +
                               treeName + ";" + std::to_string(cyc);
            chain.Add(spec.c_str());
        }

        // Копируем все события в новое дерево
        outputFile->cd();
        TTree* merged = chain.CloneTree(-1, "fast");
        merged->SetName(treeName.c_str());
        merged->Write();

        log("Merged tree written to output with " +
                std::to_string(merged->GetEntries()) + " entries.",
            LogLevel::Info);

        return true;
    }
};
