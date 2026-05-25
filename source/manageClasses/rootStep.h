#pragma once

#include "TFile.h"
#include "baseStep.h"

#include <map>
#include <string>

/**
 * @brief Абстрактный класс для шагов, работающих с ROOT-файлами.
 * 
 * Открывает *все* ROOT-файлы, указанные в inputFileNames (по ролям).
 * Предоставляет доступ к каждому через map `inputFiles`.
 * Также открывает выходной файл, если указан.
 * 
 * Закрывает все открытые файлы в `finalize()`.
 */
class RootStep : public BaseStep {
public:
    RootStep(const std::string& stepName,
             const std::map<std::string, std::string>& inputFileNames,
             const std::string& outputFileName = "")
        : BaseStep(stepName, inputFileNames, outputFileName) {}

    virtual ~RootStep() = default;

protected:
    std::map<std::string, TFile*> inputFiles;  ///< Все входные ROOT-файлы по ролям
    TFile* inputFile = nullptr;                ///< Устаревшее имя: синоним для inputFiles["main"]
    TFile* outputFile = nullptr;               ///< Выходной ROOT-файл (если указан)

    virtual bool initialize() override {
        // Открываем все входные файлы
        for (const auto& [role, filename] : inputFileNames) {
            log("Opening ROOT input file [" + role + "]: " + filename, LogLevel::Debug);
            TFile* file = TFile::Open(filename.c_str(), "READ");
            if (!file || file->IsZombie()) {
                log("Failed to open input file: " + filename + " (role: " + role + ")", LogLevel::Error);
                return false;
            }
            inputFiles[role] = file;
        }

        // Синоним для старого кода
        if (inputFiles.count("main")) {
            inputFile = inputFiles["main"];
            log("ROOT file opened successfully: " + getInputFileName("main"), LogLevel::Info);
        } else {
            log("Warning: input file with role 'main' not found.", LogLevel::Warning);
        }

        // Открываем выходной файл, если указан
        if (!outputFileName.empty()) {
            outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
            if (!outputFile || outputFile->IsZombie()) {
                log("Failed to create output file: " + outputFileName, LogLevel::Error);
                return false;
            }
            log("ROOT output file created successfully: " + outputFileName, LogLevel::Info);
        }

        return true;
    }

    virtual bool finalize() override {
        // Закрываем все входные файлы
        for (auto& [role, file] : inputFiles) {
            if (file && file->IsOpen()) {
                file->Close();
                log("Closed input file [" + role + "]", LogLevel::Debug);
            }
        }

        inputFiles.clear();
        inputFile = nullptr;

        // Закрываем выходной файл
        if (outputFile && outputFile->IsOpen()) {
            outputFile->Close();
            log("Closed output file", LogLevel::Debug);
        }

        outputFile = nullptr;
        return true;
    }

    // Утилита для доступа
    TFile* getInputFile(const std::string& role) const {
        auto it = inputFiles.find(role);
        return (it != inputFiles.end()) ? it->second : nullptr;
    }
};
