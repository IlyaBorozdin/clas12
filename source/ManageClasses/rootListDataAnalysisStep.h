#pragma once

#include <memory>
#include "TTree.h"
#include "TKey.h"
#include "rootStep.h"

/**
 * @brief Базовый абстрактный класс для шагов анализа, обрабатывающих несколько TTree в ROOT-файле.
 *
 * - Загружает ROOT-файл (через RootStep)
 * - Вызывает пользовательский метод `processTree()` для каждого TTree, имя которого задаёт `nameGenerator`
 * - Используется для шагов, где логика анализа применяется ко *всем* деревьям, а не событиям внутри одного дерева
 *
 * Типичный сценарий: обрабатывать набор деревьев, полученных после разбиения по секторам, углам и т.д.
 */
class RootListDataAnalysisStep : public RootStep {
public:
    /**
     * @param stepName Имя шага (отображается в логах)
     * @param inputFileNames Словарь входных файлов (ожидается, что файл с ролью "main" содержит TTree)
     * @param nameGenerator Функция, выдающая последовательность имён деревьев для анализа
     * @param outputFileName Имя выходного файла (если нужно)
     */
    RootListDataAnalysisStep(const std::string& stepName,
                             const std::map<std::string, std::string>& inputFileNames,
                             std::function<bool(std::string& nextName)> nameGenerator,
                             const std::string& outputFileName = "")
        : RootStep(stepName, inputFileNames, outputFileName),
          nameGenerator(nameGenerator) {}

    virtual ~RootListDataAnalysisStep() = default;

protected:
    std::function<bool(std::string& nextName)> nameGenerator;
    size_t totalTrees = 0;     ///< Общее количество деревьев
    size_t treesSkipped = 0;

    /**
     * @brief Инициализация: открывает файл и подсчитывает количество TTree.
     * 
     * Метод полезен в логах и для отладки. Само чтение деревьев происходит в run().
     */
    virtual bool initialize() override {
        if (!RootStep::initialize()) {
            return false;
        }

        totalTrees = 0;
        TIter next(inputFile->GetListOfKeys());
        TKey* key;

        while ((key = (TKey*)next())) {
            if (std::string(key->GetClassName()) == "TTree") {
                ++totalTrees;
            }
        }

        log("Discovered " + std::to_string(totalTrees) + " TTree object(s) in the input file.", LogLevel::Info);

        if (totalTrees == 0) {
            log("Warning: no TTree objects found in the input file.", LogLevel::Warning);
        }

        return true;
    }


    /**
     * @brief Основной метод выполнения: вызывает обработку каждого дерева по имени.
     * 
     * Использует генератор имён деревьев (nameGenerator) и вызывает для каждого:
     * - метод processTree()
     */
    virtual bool run() override {
        if (!inputFile) {
            log("Input file is not open", LogLevel::Error);
            return false;
        }

        std::string name;

        while (nameGenerator(name)) {
            TTree* tree = nullptr; 
            inputFile->GetObject(name.c_str(), tree);

            if (!tree) {
                ++treesSkipped;
                continue;
            }

            logProgress();
            processTree(tree);
        }
        log("Finished processing all requested trees.", LogLevel::Info);

        return true;
    }

    virtual bool finalize() override {
        if (totalTrees == 0) {
            log("No trees were requested for processing.", LogLevel::Warning);
        } else {
            log("Final stats: " + std::to_string(totalTrees - treesSkipped) + " tree(s) processed, "
                + std::to_string(treesSkipped) + " skipped (" +
                std::to_string(100.0 * treesSkipped / totalTrees) + "%).", LogLevel::Info);
        }

        return true;
    }


    /**
     * @brief Метод обработки одного дерева. Вызывается один раз на дерево.
     * 
     * Реализуется в производном классе. Содержит логику анализа одного TTree.
     */
    virtual void processTree(TTree* tree) = 0;

    /**
     * @brief Метод логирования прогресса.
     * 
     * Может вызываться на каждой итерации `run()` для отображения прогресса.
     * Реализуется в потомке (например, логгирование каждого N-го дерева).
     */
    virtual void logProgress() = 0;
};
