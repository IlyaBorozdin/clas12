#pragma once

#include "source/ManageClasses/baseStep.h" 
#include "choiceRoot.h"

/**
 * @brief Шаг конвертации списка HIPO-файлов в ROOT-формат с использованием класса ChoiceRoot.
 * 
 * Данный шаг оборачивает существующий класс `ChoiceRoot`, чтобы встроить его
 * в модульную архитектуру на базе `BaseStep`, сохранив существующую логику анализа.
 * 
 * Входной файл — текстовый список путей к HIPO-файлам ("hipo_list").
 * Выходной файл — ROOT-файл с результатами (`TTree`), создаваемый ChoiceRoot.
 */
class HipoConversionStep : public BaseStep {
public:
    /**
     * @param stepName имя шага (для логирования)
     * @param inputFileName имя входного TXT-файла
     * @param outputFileName имя выходного ROOT-файла
     */
    HipoConversionStep(const std::string& stepName,
                          const std::string& inputFileName,
                          const std::string& outputFileName,
                          bool appendMode = false,
                          const std::string& treeName = "ExpData")
        : BaseStep(stepName, { {"hipo_list", inputFileName} }, outputFileName, appendMode),
          choice(nullptr),
          appendMode(appendMode),
          treeName(treeName) {}

    ~HipoConversionStep() override {
        delete choice;
    }

    void describe() const override {
        BaseStep::describe();
        log("This step wraps ChoiceRoot to convert HIPO data into ROOT tree.", LogLevel::Debug);
    }

protected:
    bool initialize() override {
        log("Initializing HipoConversionStep...", LogLevel::Info);

        if (!hasInputFile("hipo_list")) {
            log("Missing input file with role 'hipo_list'", LogLevel::Error);
            return false;
        }

        const std::string& hipoFileList = getInputFileName("hipo_list");
        const std::string& output = getOutputFileName();

        log("Creating ChoiceRoot with input: " + hipoFileList + ", output: " + output, LogLevel::Debug);
        choice = new ChoiceRoot(hipoFileList.c_str(), output.c_str(), appendMode, treeName.c_str());

        return true;
    }

    bool run() override {
        log("Running ChoiceRoot::analysisCycle()...", LogLevel::Info);
        choice->analysisCycle();
        log("ChoiceRoot finished successfully.", LogLevel::Info);
        return true;
    }

    bool finalize() override {
        log("Finalizing HipoConversionStep.", LogLevel::Info);
        return true;
    }

private:
    ChoiceRoot* choice;

    bool appendMode;
    std::string treeName;
};
