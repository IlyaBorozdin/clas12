#pragma once

#include "analysisBranch.h"

/**
 * @brief Контроллер всего анализа. Управляет набором веток анализа.
 */
class AnalysisManager {
public:
    void addBranch(const std::string& name, const std::shared_ptr<AnalysisBranch>& branch) {
        branches[name] = branch;
    }

    bool runBranch(const std::string& name) {
        auto it = branches.find(name);
        if (it != branches.end()) {
            std::cout << "[MANAGER] Running branch: " << name << std::endl;
            bool success = it->second->runAll();
            if (!success) {
                std::cerr << "[MANAGER] Branch " << name << " finished with errors." << std::endl;
            }
            return success;
        } else {
            std::cerr << "[MANAGER] Branch not found: " << name << std::endl;
            return false;
        }
    }

    void describe() const {
        std::cout << "[MANAGER] Contains " << branches.size() << " branches:" << std::endl;
        for (const auto& [name, branch] : branches) {
            std::cout << "\nBranch: " << name << std::endl;
            branch->describe();
        }
    }

    bool hasBranch(const std::string& name) const {
        return branches.find(name) != branches.end();
    }

    void removeBranch(const std::string& name) {
        branches.erase(name);
    }

private:
    std::map<std::string, std::shared_ptr<AnalysisBranch>> branches;
};
