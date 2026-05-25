#pragma once

#include "source/manageClasses/rootTreeStep.h"
#include "utils/fitModel.h"

class GenericFitAnalysisStep : public RootTreeStep {
protected:
    std::unique_ptr<FitModel> model;
    
    // Общие кинематические переменные (остаются здесь)
    double yield, eYield, chi2ndf, width_mm, downEdge, upEdge;
    double neutronIntegral, neutronValue, neutronPosition, maxValue, maxPosition;
    int index_q2, index_w, index_cos_theta, index_phi;

public:
    GenericFitAnalysisStep(const std::string& stepName,
                           const std::map<std::string, std::string>& inputFileNames,
                           const std::string& outputFileName,
                           std::unique_ptr<FitModel> fitModel)
        : RootTreeStep(stepName, inputFileNames, "FitData", outputFileName),
          model(std::move(fitModel)) {}

    void setBranchAddresses() override {
        // 1. Делегируем модели привязку её специфичных параметров
        model->setupBranches(inputTree);

        // 2. Сами привязываем общие переменные
        inputTree->SetBranchAddress("yield", &yield);
        inputTree->SetBranchAddress("eYield", &eYield);
        inputTree->SetBranchAddress("chi2ndf", &chi2ndf);

        inputTree->SetBranchAddress("downEdge", &downEdge);
        inputTree->SetBranchAddress("upEdge", &upEdge);
        inputTree->SetBranchAddress("width_mm", &width_mm);
        inputTree->SetBranchAddress("neutronIntegral", &neutronIntegral);
        inputTree->SetBranchAddress("neutronValue", &neutronValue);
        inputTree->SetBranchAddress("neutronPosition", &neutronPosition);
        inputTree->SetBranchAddress("maxValue", &maxValue);
        inputTree->SetBranchAddress("maxPosition", &maxPosition);

        inputTree->SetBranchAddress("index_q2", &index_q2);
        inputTree->SetBranchAddress("index_w", &index_w);
        inputTree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        inputTree->SetBranchAddress("index_phi", &index_phi);
    }

    // run() теперь тоже универсален
    bool run() override {
        const Long64_t nEntries = inputTree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            inputTree->GetEntry(i);
            if (check()) {
                processEvent();
            }
        }
        return true;
    }

    virtual bool check() { return true; }

    virtual void processEvent() = 0;
};