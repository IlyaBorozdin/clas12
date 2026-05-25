#pragma once

#include <memory>
#include <vector>
#include <algorithm>

#include "baseStep.h"

/**
 * @brief Координатор шагов в рамках одной ветки анализа
 *        Управляет последовательностью шагов (BaseStep), вызывая их по порядку
 */
class AnalysisBranch {
public:
    void addStep(const std::shared_ptr<BaseStep>& step) {
        steps.push_back(step);
    }

    bool runAll() {
        for (auto& step : steps) {
            step->safeRun();
        }
        return allDone();
    }

    void describe() const {
        std::cout << "[BRANCH] Contains " << steps.size() << " steps:" << std::endl;
        for (const auto& step : steps) {
            std::cout << "  - " << step->getStepName() << " [Status: " << to_string(step->getStatus()) << "]\n";
        }
    }

    const std::vector<std::shared_ptr<BaseStep>>& getSteps() const {
        return steps;
    }

    std::shared_ptr<BaseStep> getStepByName(const std::string& name) const {
        for (const auto& step : steps) {
            if (step->getStepName() == name) return step;
        }
        return nullptr;
    }

    bool allDone() const {
        return std::all_of(steps.begin(), steps.end(), [](const auto& step) {
            return step->getStatus() == StepStatus::Done || step->getStatus() == StepStatus::Skipped;
        });
    }

    // Итераторы для внешней итерации
    auto begin() const { return steps.begin(); }
    auto end() const { return steps.end(); }

private:
    std::vector<std::shared_ptr<BaseStep>> steps;
};
