#pragma once

#include <vector>

#include "TTree.h"
#include "TKey.h"

#include "rootStep.h"

/**
 * @brief Абстрактный класс для шагов, обрабатывающих несколько TTree в одном ROOT-файле.
 * Наследуется от RootStep. Загружает все деревья из файла и предоставляет доступ к ним.
 */
class RootListTreesStep : public RootStep {
public:
    RootListTreesStep(const std::string& stepName,
                      const std::map<std::string, std::string>& inputFileNames,
                      const std::string& outputFileName = "")
        : RootStep(stepName, inputFileNames, outputFileName) {}

    ~RootListTreesStep() override {
        for (TTree* tree : trees) {
            delete tree;
        }
    }

protected:
    std::vector<TTree*> trees;  ///< Все загруженные деревья
    size_t totalTrees = 0;

    /**
     * @brief Метод инициализации: находит все TTree в файле.
     */
    bool initialize() override {
        if (!RootStep::initialize()) {
            return false;
        }

        TIter next(inputFile->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TTree* tree = dynamic_cast<TTree*>(key->ReadObj());
            if (!tree) {
                log("Skipping non-TTree object: " + std::string(key->GetName()), LogLevel::Warning);
                continue;
            }

            tree->SetDirectory(nullptr); // Отвязываем от файла
            trees.push_back(tree);
        }

        totalTrees = trees.size();
        log("Found " + std::to_string(totalTrees) + " TTree objects in input file.", LogLevel::Info);

        if (totalTrees == 0) {
            log("No TTrees found in the input file.", LogLevel::Warning);
        }

        return true;
    }

    /**
     * @brief Возвращает дерево по индексу (если в пределах диапазона), иначе nullptr
     */
    TTree* getTreeByIndex(size_t index) const {
        return (index < trees.size()) ? trees[index] : nullptr;
    }

    TTree* findTreeByName(const std::string& name) const {
        for (TTree* tree : trees) {
            if (name == tree->GetName()) return tree;
        }
        return nullptr;
    }
};
