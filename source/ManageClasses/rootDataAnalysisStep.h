#pragma once

#include "rootTreeStep.h"

/**
 * @brief Базовый класс для шагов анализа данных на основе событий в ROOT-дереве.
 *
 * Наследует RootTreeStep и реализует типовой цикл обработки:
 * - перебор всех записей дерева;
 * - вызов `processEvent()` для каждой записи;
 * - ведение логов с помощью `logProgress()`.
 *
 * Подходит для шагов анализа событий с известной фиксированной структурой (например: e_*, pi_*).
 * Ветка `setBranchAddresses()` связывает стандартные переменные с ветками дерева.
 */
class RootDataAnalysisStep : public RootTreeStep {
public:
    /**
     * @param stepName       Имя шага (используется в логах)
     * @param inputFileName  Имя входного ROOT-файла (предполагается, что содержит дерево "ExpData")
     * @param outputFileName Имя выходного файла (если есть)
     */
    RootDataAnalysisStep(const std::string& stepName,
                         const std::string& inputFileName,
                         const std::string& outputFileName = "")
        : RootTreeStep(stepName, { {"main", inputFileName} }, "ExpData", outputFileName) {}

    virtual ~RootDataAnalysisStep() = default;

protected:
    // --- Переменные, связанные с ветками входного дерева (Input TTree) ---
    double e_px = 0, e_py = 0, e_pz = 0;     ///< Компоненты импульса электрона
    double pi_px = 0, pi_py = 0, pi_pz = 0;  ///< Компоненты импульса пиона


    /**
     * @brief Реализация основного алгоритма: цикл по всем событиям дерева.
     * 
     * - вызывает initialize() (открытие файла, получение дерева, связывание веток);
     * - проходит по всем записям дерева, вызывая processEvent() для каждой;
     * - логирует прогресс;
     * - завершает шаг с finalize().
     * 
     * @return true, если шаг выполнен успешно, иначе false
     */
    virtual bool run() override {
        const Long64_t nEntries = inputTree->GetEntries();
        log("Starting event loop: " + std::to_string(nEntries) + " entries.", LogLevel::Info);

        for (Long64_t i = 0; i < nEntries; ++i) {
            inputTree->GetEntry(i);
            logProgress();
            processEvent();
        }
        log("Event loop completed.", LogLevel::Info);

        return true;
    }

    /**
     * @brief Привязывает переменные-члены к соответствующим веткам в дереве.
     * 
     * Данный метод реализует привязку к фиксированному набору веток:
     * - импульс электрона (e_px, e_py, e_pz)
     * - импульс пиона (pi_px, pi_py, pi_pz)
     * 
     * При необходимости анализа других структур — переопределите метод.
     */
    virtual void setBranchAddresses() override {
        inputTree->SetBranchAddress("e_px", &e_px);
        inputTree->SetBranchAddress("e_py", &e_py);
        inputTree->SetBranchAddress("e_pz", &e_pz);
        inputTree->SetBranchAddress("pi_px", &pi_px);
        inputTree->SetBranchAddress("pi_py", &pi_py);
        inputTree->SetBranchAddress("pi_pz", &pi_pz);
    }

    /**
     * @brief Метод обработки одного события.
     * 
     * Реализуется в наследнике. Использует переменные, связанные с ветками дерева.
     */
    virtual void processEvent() = 0;

    /**
     * @brief Логгирует прогресс обработки событий (например, каждые N событий).
     * 
     * Реализуется в наследнике. Может использовать totalEvents или текущий счётчик.
     */
    virtual void logProgress() = 0;

};
