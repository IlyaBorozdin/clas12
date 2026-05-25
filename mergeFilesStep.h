#pragma once

#include "source/manageClasses/rootStep.h"

#include "TChain.h"
#include "TKey.h"
#include "TTree.h"

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

/**
 * @brief Шаг для объединения последних (current) циклов дерева из нескольких входных ROOT-файлов.
 *
 * Использует уже открытые файлы из RootStep::inputFiles.
 * Для каждого входного файла берёт последний цикл дерева `treeName` и добавляет его в TChain.
 * Затем копирует все события в единое дерево в outputFile.
 */
class MergeFilesStep : public RootStep {
public:
    MergeFilesStep(const std::string& stepName,
                   const std::map<std::string, std::string>& inputFilesMap,
                   const std::string& outputFileName,
                   const std::string& treeName = "ExpData")
        : RootStep(stepName, inputFilesMap, outputFileName),
          treeName(treeName) {}

protected:
    std::string treeName;

    bool run() override {
        if (inputFiles.empty()) {
            log("No input files available (inputFiles is empty). Did RootStep::initialize succeed?", LogLevel::Error);
            return false;
        }

        // Собираем TChain из последних циклов treeName в каждом открытом файле
        TChain chain(treeName.c_str());
        Long64_t expectedTotal = 0;
        int addedFiles = 0;

        for (const auto& [role, filePtr] : inputFiles) {
            TFile* file = filePtr;
            if (!file) {
                log("Null TFile* for role: " + role, LogLevel::Warning);
                continue;
            }
            if (file->IsZombie()) {
                log("TFile is a Zombie for role: " + role + " (skipping)", LogLevel::Warning);
                continue;
            }

            TList* keys = file->GetListOfKeys();
            if (!keys || keys->IsEmpty()) {
                log("No keys in file for role: " + role + " (file: " + std::string(file->GetName()) + ")", LogLevel::Warning);
                continue;
            }

            // Найдём все циклы для дерева treeName
            std::vector<int> cycles;
            for (TObject* obj : *keys) {
                auto* key = dynamic_cast<TKey*>(obj);
                if (!key) continue;
                if (treeName == std::string(key->GetName())) {
                    cycles.push_back(key->GetCycle());
                }
            }

            if (cycles.empty()) {
                log("No cycles for tree '" + treeName + "' in file: " + std::string(file->GetName()), LogLevel::Warning);
                continue;
            }

            int maxCycle = *std::max_element(cycles.begin(), cycles.end());

            // Формируем спецификацию для TChain: "filename/TreeName;cycle"
            std::string fname = file->GetName(); // имя файла на диске (как открытое)
            std::string spec = fname + "/" + treeName + ";" + std::to_string(maxCycle);
            chain.Add(spec.c_str());
            ++addedFiles;

            // Попробуем получить указатель на конкретный цикл дерева чтобы узнать количество записей
            std::string keyName = treeName + ";" + std::to_string(maxCycle);
            TTree* tmp = static_cast<TTree*>(file->Get(keyName.c_str()));
            Long64_t entries = tmp ? tmp->GetEntries() : 0;
            expectedTotal += entries;

            log("Added " + keyName + " from " + fname + " (" + std::to_string(entries) + " entries).", LogLevel::Info);
        }

        if (chain.GetNtrees() == 0) {
            log("No trees were added to the chain (no valid last-cycle trees found).", LogLevel::Error);
            return false;
        }

        // Клонируем все события в единый TTree в outputFile
        outputFile->cd();
        TTree* merged = chain.CloneTree(-1, "fast");
        if (!merged) {
            log("CloneTree returned null.", LogLevel::Error);
            return false;
        }
        merged->SetName(treeName.c_str());
        merged->Write();

        log("Merged tree written to output with " + std::to_string(merged->GetEntries()) +
            " entries (expected total from per-file checks: " + std::to_string(expectedTotal) + ").", LogLevel::Info);

        log("Files processed: " + std::to_string(addedFiles) + ", chain ntrees: " + std::to_string(chain.GetNtrees()), LogLevel::Debug);

        return true;
    }
};
