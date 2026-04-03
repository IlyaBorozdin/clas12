#pragma once

#include "source/ManageClasses/baseStep.h"
#include "hipoStatsBase.h"

/**
 * @brief Шаг сбора статистической информации из HIPO-файлов.
 * * Оборачивает классы, наследуемые от HipoStatsBase, позволяя интегрировать
 * сбор статистики в общий пайплайн BaseStep.
 */
template <typename StatsType>
class HipoStatsStepT : public BaseStep {
public:
    /**
     * @param stepName имя шага
     * @param inputFileName текстовый список HIPO-файлов
     * @param outputFileName текстовый файл отчета (.txt)
     */
    HipoStatsStepT(const std::string& stepName,
                   const std::string& inputFileName,
                   const std::string& outputFileName)
        : BaseStep(stepName, { {"hipo_list", inputFileName} }, outputFileName),
          statsProcessor(nullptr) {}

    ~HipoStatsStepT() override {
        delete statsProcessor;
    }

    void describe() const override {
        BaseStep::describe();
        log("This step wraps a HipoStatsBase derivative to generate statistical reports.", LogLevel::Debug);
    }

protected:
    bool initialize() override {
        log("Initializing HipoStatsStep...", LogLevel::Info);

        if (!hasInputFile("hipo_list")) {
            log("Missing input file with role 'hipo_list'", LogLevel::Error);
            return false;
        }

        const std::string& hipoFileList = getInputFileName("hipo_list");
        const std::string& outputReport = getOutputFileName();

        log("Creating StatsProcessor with input: " + hipoFileList + ", report: " + outputReport, LogLevel::Debug);
        
        // Создаем объект статистики (например, SimpleParticleStats)
        statsProcessor = new StatsType(hipoFileList.c_str(), outputReport.c_str());

        return (statsProcessor != nullptr);
    }

    bool run() override {
        log("Running statistics analysis cycle...", LogLevel::Info);
        
        // Вызываем метод базового класса HipoDataAnalysis через нашего наследника
        statsProcessor->analysisCycle();
        
        log("Statistics analysis finished successfully.", LogLevel::Info);
        return true;
    }

    bool finalize() override {
        log("Finalizing HipoStatsStep.", LogLevel::Info);
        return true;
    }

private:
    StatsType* statsProcessor;
};

using SimpleParticleStatsStep = HipoStatsStepT<SimpleParticleStats>;