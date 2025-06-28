#pragma once

#include "TTree.h"

#include "rootStep.h"

/**
 * @brief Абстрактный класс для шагов, работающих с одним TTree в ROOT-файле.
 * 
 * Наследует RootStep и реализует:
 * - автоматическое открытие ROOT-файла по роли "main"
 * - извлечение дерева по имени (по умолчанию "ExpData")
 * - получение количества событий
 * - вызов абстрактного метода setBranchAddresses()
 * 
 * Подходит для шагов, обрабатывающих структуру событий в одном дереве.
 */
class RootTreeStep : public RootStep {
public:
    /**
     * @param stepName Имя шага (для логирования и отладки)
     * @param inputFileNames Карта входных файлов, должен быть задан файл с ролью "main"
     * @param inputTreeName Имя дерева в ROOT-файле (по умолчанию — "ExpData")
     * @param outputFileName Имя выходного файла (опционально)
     */
    RootTreeStep(const std::string& stepName,
                 const std::map<std::string, std::string>& inputFileNames,
                 const std::string& inputTreeName = "ExpData",
                 const std::string& outputFileName = "")
        : RootStep(stepName, inputFileNames, outputFileName),
          inputTreeName(inputTreeName) {}

    /// Деструктор: закрывает и удаляет входной файл, если он был открыт
    virtual ~RootTreeStep() = default;

protected:
    TTree* inputTree = nullptr;          ///< Указатель на дерево, извлечённое из ROOT-файла (inputFile)
    std::string inputTreeName;           ///< Имя дерева, которое будет загружено (по умолчанию "ExpData")
    size_t totalEvents = 0;              ///< Общее число событий в дереве (inputTree->GetEntries())


    /**
     * @brief Метод инициализации: находит дерево по имени, настраивает ветки.
     * 
     * @return true если инициализация прошла успешно, иначе false
     */
    virtual bool initialize() override {
        if (!RootStep::initialize()) {
            return false;
        }

        log("Looking for TTree: " + inputTreeName, LogLevel::Debug);

        inputFile->GetObject(inputTreeName.c_str(), inputTree);
        if (!inputTree) {
            log("TTree '" + inputTreeName + "' not found in input file.", LogLevel::Error);
            return false;
        }

        totalEvents = inputTree->GetEntries();
        if (totalEvents == 0) {
            log("TTree '" + inputTreeName + "' is empty.", LogLevel::Warning);
        } else {
            log("TTree '" + inputTreeName + "' loaded. Total events: " + std::to_string(totalEvents), LogLevel::Info);
        }

        setBranchAddresses();
        log("Branch addresses set successfully.", LogLevel::Debug);
        return true;
    }


    /**
     * @brief Абстрактный метод: связывает переменные с ветками TTree.
     * 
     * Каждый наследник должен реализовать этот метод в соответствии со структурой дерева.
     * Обычно включает вызовы вида inputTree->SetBranchAddress("var", &varPtr).
     */
    virtual void setBranchAddresses() = 0;
};
