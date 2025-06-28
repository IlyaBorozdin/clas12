#pragma once

#include "TFile.h"

#include "baseStep.h"

/**
 * @brief Абстрактный класс для шагов, работающих с ROOT-файлом.
 * 
 * Открывает ROOT-файл, указанный как "main" во входных файлах.
 * Предоставляет доступ к открытому TFile* в производных классах.
 * Закрывает файл в деструкторе.
 * 
 * Используется как базовый класс для шагов, работающих с деревьями, гистограммами и т.д.
 */
class RootStep : public BaseStep {
public:
    /**
     * @param stepName Имя шага (для логирования и отладки)
     * @param inputFileNames Карта входных файлов (ожидается наличие файла с ролью "main")
     * @param outputFileName Имя выходного файла (опционально)
     */
    RootStep(const std::string& stepName,
             const std::map<std::string, std::string>& inputFileNames,
             const std::string& outputFileName = "")
        : BaseStep(stepName, inputFileNames, outputFileName) {}

    /// @brief Деструктор. Закрывает и освобождает входной ROOT-файл, если он был открыт.
    /// Это освобождает все связанные ресурсы и предотвращает утечки памяти.
    virtual ~RootStep() {
        if (inputFile && inputFile->IsOpen()) {
            inputFile->Close();
            delete inputFile;
        }
    }

protected:
    TFile* inputFile = nullptr; ///< Указатель на входной ROOT-файл. Открывается в initialize(), закрывается в деструкторе.

    /**
     * @brief Метод инициализации: открывает входной файл.
     * 
     * @return true если инициализация прошла успешно, иначе false
     */
    virtual bool initialize() override {
        // Проверка наличия входного файла с ролью "main"
        if (!hasInputFile("main")) {
            log("No input file with role 'main'.", LogLevel::Error);
            return false;
        }

        // Получение имени файла
        const std::string& fileName = getInputFileName("main");

        // Открытие ROOT-файла на чтение
        log("Opening ROOT file: " + fileName, LogLevel::Debug);
        inputFile = TFile::Open(fileName.c_str(), "READ");
        if (!inputFile || inputFile->IsZombie()) {
            log("Failed to open file: " + fileName, LogLevel::Error);
            return false;
        }
        log("ROOT file opened successfully: " + fileName, LogLevel::Info);

        return true;
    }
};
