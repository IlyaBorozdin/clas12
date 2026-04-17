#pragma once

#include <memory>
#include <functional>
#include "TKey.h"
#include "TFile.h"
#include "TClass.h"

#include "rootStep.h"

/**
 * @brief Абстрактный шаблонный класс для шагов анализа, обрабатывающих объекты одного типа из ROOT-файла.
 *
 * Параметризован типом ROOT-объекта (например, TTree или TH1F).
 */
template <typename T>
class RootListObjectAnalysisStep : public RootStep {
public:
    RootListObjectAnalysisStep(const std::string& stepName,
                               const std::map<std::string, std::string>& inputFileNames,
                               std::function<bool(std::string&)> nameGenerator,
                               const std::string& outputFileName = "")
        : RootStep(stepName, inputFileNames, outputFileName),
          nameGenerator(nameGenerator) {}

    virtual ~RootListObjectAnalysisStep() = default;

protected:
    std::function<bool(std::string&)> nameGenerator;
    size_t totalObjects = 0;
    size_t objectsSkipped = 0;

    virtual bool initialize() override {
        if (!RootStep::initialize())
            return false;

        // Подсчитаем количество объектов нужного типа
        TIter next(inputFile->GetListOfKeys());
        TKey* key;
        TClass* expectedClass = T::Class();

        while ((key = (TKey*)next())) {
            if (std::string(key->GetClassName()) == expectedClass->GetName()) {
                ++totalObjects;
            }
        }

        log("Discovered " + std::to_string(totalObjects) + " object(s) of type " + std::string(expectedClass->GetName()), LogLevel::Info);

        if (totalObjects == 0) {
            log("Warning: no objects of type " + std::string(expectedClass->GetName()) + " found in the input file.", LogLevel::Warning);
        }

        return true;
    }

    virtual bool run() override {
        if (!inputFile) {
            log("Input file is not open", LogLevel::Error);
            return false;
        }

        std::string name;
        while (nameGenerator(name)) {
            T* object = nullptr;
            inputFile->GetObject(name.c_str(), object);

            if (!object) {
                ++objectsSkipped;
                continue;
            }

            logProgress();
            if (check(object)) {
                processObject(object);
            }
        }

        log("Finished processing all requested objects.", LogLevel::Info);
        return true;
    }

    virtual bool finalize() override {
        if (totalObjects == 0) {
            log("No objects were requested for processing.", LogLevel::Warning);
        } else {
            log("Final stats: " + std::to_string(totalObjects - objectsSkipped) + " processed, "
                + std::to_string(objectsSkipped) + " skipped (" +
                std::to_string(100.0 * objectsSkipped / totalObjects) + "%).", LogLevel::Info);
        }
        return RootStep::finalize();
    }

    virtual bool check(T* object) { return true; }

    // 👇 Специализируется в потомке
    virtual void processObject(T* object) = 0;
    virtual void logProgress() {};
};
