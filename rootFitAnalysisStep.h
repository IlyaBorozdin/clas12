#pragma once

#include "source/ManageClasses/rootTreeStep.h"

class RootFitAnalysisStep : public RootTreeStep {
public:
    RootFitAnalysisStep(const std::string& stepName,
                        const std::map<std::string, std::string>& inputFileNames,
                        const std::string& outputFileName = "")
        : RootTreeStep(stepName, inputFileNames, "FitData", outputFileName) {}

    virtual ~RootFitAnalysisStep() = default;

protected:
    double ampGauss, meanGauss, stdDevGauss, asymGauss, ampGaussDelta, meanGaussDelta, stdDevGaussDelta, paramA, paramB, paramC, width_mm;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss, eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta, errA, errB, errC;
    double downEdge, upEdge;
    int index_q2, index_w, index_cos_theta, index_phi;

    virtual bool run() override {
        const Long64_t nEntries = inputTree->GetEntries();
        log("Starting event loop: " + std::to_string(nEntries) + " entries.", LogLevel::Info);

        for (Long64_t i = 0; i < nEntries; ++i) {
            inputTree->GetEntry(i);
            logProgress();
            if (check()) {
                processEvent();
            }
        }
        log("Event loop completed.", LogLevel::Info);

        return true;
    }

    void setBranchAddresses() override {
        inputTree->SetBranchAddress("ampGauss", &ampGauss);
        inputTree->SetBranchAddress("meanGauss", &meanGauss);
        inputTree->SetBranchAddress("stdDevGauss", &stdDevGauss);
        inputTree->SetBranchAddress("asymGauss", &asymGauss);
        inputTree->SetBranchAddress("ampGaussDelta", &ampGaussDelta);
        inputTree->SetBranchAddress("meanGaussDelta", &meanGaussDelta);
        inputTree->SetBranchAddress("stdDevGaussDelta", &stdDevGaussDelta);
        inputTree->SetBranchAddress("paramA", &paramA);
        inputTree->SetBranchAddress("paramB", &paramB);
        inputTree->SetBranchAddress("paramC", &paramC);
        inputTree->SetBranchAddress("eAmpGauss", &eAmpGauss);
        inputTree->SetBranchAddress("eMeanGauss", &eMeanGauss);
        inputTree->SetBranchAddress("eStdDevGauss", &eStdDevGauss);
        inputTree->SetBranchAddress("eAsymGauss", &eAsymGauss);
        inputTree->SetBranchAddress("eAmpGaussDelta", &eAmpGaussDelta);
        inputTree->SetBranchAddress("eMeanGaussDelta", &eMeanGaussDelta);
        inputTree->SetBranchAddress("eStdDevGaussDelta", &eStdDevGaussDelta);
        inputTree->SetBranchAddress("errA", &errA);
        inputTree->SetBranchAddress("errB", &errB);
        inputTree->SetBranchAddress("errC", &errC);

        inputTree->SetBranchAddress("downEdge", &downEdge);
        inputTree->SetBranchAddress("upEdge", &upEdge);
        inputTree->SetBranchAddress("width_mm", &width_mm);
        inputTree->SetBranchAddress("index_q2", &index_q2);
        inputTree->SetBranchAddress("index_w", &index_w);
        inputTree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        inputTree->SetBranchAddress("index_phi", &index_phi);
    }

    virtual bool check() { return true; }

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