#pragma once

#include <string>
#include <map>
#include <iostream>
#include <fstream>

// Перечисление статусов выполнения шага
enum class StepStatus { NotStarted, Running, Skipped, Done, Failed };

// Уровни логирования: задают важность сообщения и влияют на формат вывода.
// Можно использовать для фильтрации логов (например, выводить только Warning+Error).
enum class LogLevel { Info, Warning, Error, Debug };

inline std::string to_string(StepStatus status) {
    switch (status) {
        case StepStatus::NotStarted: return "NotStarted";
        case StepStatus::Running:    return "Running";
        case StepStatus::Skipped:    return "Skipped";
        case StepStatus::Done:       return "Done";
        case StepStatus::Failed:     return "Failed";
    }
    return "Unknown";
}

inline std::string logPrefix(LogLevel level) {
    switch (level) {
        case LogLevel::Info:    return "[INFO]   ";
        case LogLevel::Warning: return "[WARNING]";
        case LogLevel::Error:   return "[ERROR]  ";
        case LogLevel::Debug:   return "[DEBUG]  ";
    }
    return "[LOG]";
}

/**
 * @brief Абстрактный базовый класс для шагов обработки данных.
 * Каждый шаг должен реализовать метод run(), а также initialize() и finalize().
 * 
 * Предоставляет механизм для:
 * - управления входными/выходными файлами
 * - отслеживания статуса выполнения
 * - безопасного запуска шага с логированием и отловом исключений
 */
class BaseStep {
public:
    /**
     * @param stepName имя шага (для логирования и идентификации)
     * @param inputFileNames словарь входных файлов: роль → имя файла
     * @param outputFileName имя выходного файла
     */
    BaseStep(const std::string& stepName,
             const std::map<std::string, std::string>& inputFileNames = {},
             const std::string& outputFileName = "",
             bool forceExecution = false)
        : stepName(stepName),
        inputFileNames(inputFileNames),
        outputFileName(outputFileName),
        forceExecution(forceExecution) {}


    virtual ~BaseStep() = default;

    /// Проверяет, что все входные файлы существуют
    virtual bool isValid() const {
        return allInputFilesExist();
    }

    /// Проверяет, существует ли выходной файл (т.е. шаг уже выполнен)
    virtual bool isDone() const {
        return !forceExecution && outputFileExists();
    }

    const std::string& getStepName() const { return stepName; }
    StepStatus getStatus() const { return status; }
    const std::string& getOutputFileName() const { return outputFileName; }

    /// Проверка наличия входного файла по его роли
    bool hasInputFile(const std::string& role) const {
        return inputFileNames.find(role) != inputFileNames.end();
    }

    /// Получение имени входного файла по роли (бросает исключение при отсутствии)
    const std::string& getInputFileName(const std::string& role) const {
        return inputFileNames.at(role);
    }

    /// Назначение входного файла по роли
    void setInputFile(const std::string& role, const std::string& filename) {
        inputFileNames[role] = filename;
    }

    /// Установка имени выходного файла
    void setOutputFile(const std::string& filename) {
        outputFileName = filename;
    }

    // Метки статуса (для управления изнутри шага)
    void markAsRunning()     { status = StepStatus::Running; }
    void markAsDone()        { status = StepStatus::Done; }
    void markAsSkipped()     { status = StepStatus::Skipped; }
    void markAsFailed()      { status = StepStatus::Failed; }
    void markAsNotStarted()  { status = StepStatus::NotStarted; }

    /// Описание шага: логирует список входных и выходного файлов, текущий статус
    virtual void describe() const {
        log("Step description:");
        if (!inputFileNames.empty()) {
            log("Input files:", LogLevel::Debug);
            for (const auto& [role, name] : inputFileNames) {
                log(" - " + role + ": " + name, LogLevel::Debug);
            }
        } else {
            log("No input files specified.", LogLevel::Warning);
        }
        log("Output file: " + outputFileName, LogLevel::Debug);
        log("Status: " + to_string(status), LogLevel::Info);
    }

    /**
     * @brief Безопасный запуск шага.
     * 
     * Включает:
     * - проверку, не был ли шаг уже выполнен (`isDone`)
     * - проверку существования всех входных файлов (`isValid`)
     * - вызов initialize(), run(), finalize()
     * - отлов исключений
     * 
     * Изменяет внутренний статус (Running → Done/Failed/Skipped).
     * Используется как стандартный способ запуска шага.
     */
    void safeRun() {
        markAsRunning();
        log("Step started.", LogLevel::Info);

        try {
            if (isDone()) {
                log("Step is already done. Skipping.", LogLevel::Info);
                markAsSkipped();
                return;
            }
            if (!isValid()) {
                log("Step is not valid.", LogLevel::Error);
                markAsNotStarted();
                return;
            }

            describe();
            if (!initialize()) {
                log("Initialize failed.", LogLevel::Error);
                markAsFailed();
                return;
            }
            if (!run()) {
                log("Run failed.", LogLevel::Error);
                markAsFailed();
                return;
            }
            if (!finalize()) {
                log("Finalize failed.", LogLevel::Error);
                markAsFailed();
                return;
            }

            markAsDone();
            log("Step completed successfully.", LogLevel::Info);

        } catch (const std::exception& e) {
            markAsFailed();
            log(std::string("Exception caught: ") + e.what(), LogLevel::Error);
        } catch (...) {
            markAsFailed();
            log("Unknown exception caught.", LogLevel::Error);
        }
    }

protected:
    std::string stepName;                                 ///< Имя шага
    std::map<std::string, std::string> inputFileNames;    ///< Входные файлы (по ролям)
    std::string outputFileName;                           ///< Выходной файл
    StepStatus status = StepStatus::NotStarted;           ///< Текущий статус

    /**
     * @brief Основной метод выполнения шага.
     * 
     * Реализуется в производных классах.
     * Должен выполнять основную работу шага (например, обработку данных).
     * Возвращает true при успешном завершении, false — при ошибке.
     */
    virtual bool run() = 0;

    /**
     * @brief Метод инициализации.
     * 
     * Реализуется в производных классах. Вызывается перед run().
     * Используется для открытия файлов, привязки веток и подготовки структуры данных.
     * Возвращает true при успехе, false — при ошибке.
     */
    virtual bool initialize() = 0;

    /**
     * @brief Метод финализации.
     * 
     * Реализуется в производных классах. Вызывается после run().
     * Используется для сохранения результатов, закрытия файлов, очистки ресурсов.
     * Возвращает true при успехе, false — при ошибке.
     */
    virtual bool finalize() = 0;


    /// Универсальный логгер с уровнями
    void log(const std::string& message, LogLevel level = LogLevel::Info) const {
        std::cout << logPrefix(level) << " [" << stepName << "] " << message << std::endl;
    }


private:
    bool forceExecution = false;

    /// Проверка существования выходного файла
    bool outputFileExists() const {
        if (outputFileName.empty()) {
            return false;
        }
        std::ifstream outfile(outputFileName);
        return outfile.good();
    }

    /// Проверка существования всех входных файлов
    bool allInputFilesExist() const {
        for (const auto& [role, fname] : inputFileNames) {
            std::ifstream infile(fname);
            if (!infile.good()) {
                log("Missing input file for role '" + role + "': " + fname, LogLevel::Error);
                return false;
            }
        }
        return true;
    }
};
