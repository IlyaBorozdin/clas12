#pragma once

#include "source/manageClasses/rootListObjectAnalysisStep.h"
#include "utils/fitterStrategy.h"
#include "MM_project_utils.h"

#include "TH1F.h"
#include "TTree.h"
#include <memory>

class ParamFitterStepBase : public RootListObjectAnalysisStep<TH1F> {
public:
    ParamFitterStepBase(const std::string& stepName,
                        const std::string& inputFileName,
                        const std::string& outputFileName,
                        std::unique_ptr<FitterStrategy> strategy) // Принимаем любой фиттер
        : RootListObjectAnalysisStep<TH1F>(stepName, { {"main", inputFileName} }, cellGenerator, outputFileName), 
          fitter(std::move(strategy)) {}

    virtual ~ParamFitterStepBase() = default;

protected:
    std::unique_ptr<FitterStrategy> fitter; // Наш универсальный движок
    TTree* outputTree = nullptr;

    // Общие переменные для всех видов фитов
    double yield, eYield, chi2ndf, width_mm, downEdge, upEdge;
    double neutronIntegral, neutronValue, neutronPosition, maxValue, maxPosition;
    int index_q2, index_w, index_cos_theta, index_phi;

    // Счётчики
    size_t processedHists = 0;
    size_t skippedHists = 0;
    size_t printEvery = 2000;

    // ЭТИ МЕТОДЫ РЕАЛИЗУЕТ НАСЛЕДНИК
    virtual void setupSpecificBranches() = 0;   // Создать ветки для конкретных параметров
    virtual void fillSpecificParams(TF1* f) = 0; // Заполнить их из TF1
    virtual void fillSpecificEmpty() = 0;       // Заполнить их дефолтами для пустых гист

    bool initialize() override {
        if (!RootListObjectAnalysisStep<TH1F>::initialize()) return false;

        outputTree = new TTree("FitData", "Tree with peaks' parameters");
        
        // Общие ветки
        outputTree->Branch("yield", &yield, "yield/D");
        outputTree->Branch("eYield", &eYield, "eYield/D");
        outputTree->Branch("chi2ndf", &chi2ndf, "chi2ndf/D");
        
        outputTree->Branch("width_mm", &width_mm, "width_mm/D");
        outputTree->Branch("downEdge", &downEdge, "downEdge/D");
        outputTree->Branch("upEdge", &upEdge, "upEdge/D");
        outputTree->Branch("neutronIntegral", &neutronIntegral, "neutronIntegral/D");
        outputTree->Branch("neutronValue", &neutronValue, "neutronValue/D");
        outputTree->Branch("neutronPosition", &neutronPosition, "neutronPosition/D");
        outputTree->Branch("maxValue", &maxValue, "maxValue/D");
        outputTree->Branch("maxPosition", &maxPosition, "maxPosition/D");

        outputTree->Branch("index_q2", &index_q2, "index_q2/I");
        outputTree->Branch("index_w", &index_w, "index_w/I");
        outputTree->Branch("index_cos_theta", &index_cos_theta, "index_cos_theta/I");
        outputTree->Branch("index_phi", &index_phi, "index_phi/I");

        setupSpecificBranches(); // Вызов специфики
        return true;
    }

    void processObject(TH1F* hist) override {
        int i, j, k, l;
        if (!parseCellName(hist->GetName(), i, j, k, l)) {
            ++skippedHists;
            return;
        }

        width_mm = hist->GetBinWidth(1);
        // downEdge = fitter->getDownEdge(hist, i, j, k, l);
        // upEdge = fitter->getUpEdge(hist, i, j, k, l);
        fitter->getBounds(hist, i, j, k, l, downEdge, upEdge);

        double events = fitter->getAnalyzer()->countEntries(hist, downEdge, upEdge);
        
        if (events < 0.5) {
            fillSpecificEmpty(); // Специфика
            yield = 0.0;
            eYield = 1.15 / width_mm;
            chi2ndf = 0.0;
            neutronPosition = 0.9396;
            neutronValue = 0.0;
            neutronIntegral = 0.0;
        } else {
            TF1* fitFunc = fitter->fit(hist, i, j, k, l);
            
            fillSpecificParams(fitFunc); // Специфика
            
            // Расчет выхода теперь делегируем фиттеру!
            yield = fitter->calculateYield(fitFunc, width_mm);
            eYield = fitter->calculateYieldError(fitFunc, width_mm);

            if (fitFunc->GetNDF() > 0) {
                chi2ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();
            } else {
                chi2ndf = 0.0; 
            }

            neutronIntegral = fitter->getAnalyzer()->getIntegral(hist, downEdge, upEdge);
            fitter->getAnalyzer()->findPeak(hist, downEdge, upEdge, neutronPosition, neutronValue);
            
            delete fitFunc;
        }

        fitter->getAnalyzer()->getGlobalMax(hist, downEdge, upEdge, maxPosition, maxValue);
        index_q2 = i; index_w = j; index_cos_theta = k; index_phi = l;

        outputTree->Fill();
    }

    bool check(TH1F* hist) override { return true; }

    bool finalize() override {
        outputFile->cd();
        outputTree->Write();
        return RootListObjectAnalysisStep<TH1F>::finalize();
    }
};