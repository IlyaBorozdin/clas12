#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"

#include "MM_project_utils.h"
#include "source/ManageClasses/rootListObjectAnalysisStep.h"
// #include "utils/histogramFitterFlexible.h"
#include "utils/histogramFitterUniversal.h"

class ParamFitterStepExt : public RootListObjectAnalysisStep<TH1F> {
public:
    ParamFitterStepExt(const std::string& stepName,
                       const std::string& inputFileName,
                       const std::string& outputFileName,
                       const std::string& configFileName = "utils/fit_configs.json")
        : RootListObjectAnalysisStep<TH1F>(stepName, { {"main", inputFileName} }, cellGenerator, outputFileName), fitter(configFileName) {}

    ~ParamFitterStepExt() override = default;

protected:
    size_t processedHists = 0;   ///< Счётчик обработанных деревьев
    size_t skippedHists = 0;     ///< Счётчик пропущенных деревьев
    size_t printEvery = 2000;    ///< Частота логирования прогресса

    bool initialize() override {
        if (!RootListObjectAnalysisStep<TH1F>::initialize()) {
            return false;
        }

        outputTree = new TTree("FitData", "Tree with peaks' parameters");
        setupBranches();

        return true;
    }

    void processObject(TH1F* hist) override {
        int i, j, k, l;
        if (!parseCellName(hist->GetName(), i, j, k, l) || fitter.getFitRegionEvents(hist, i, j, k, l) == 0) {
            ++skippedHists;
            return;
        }

        TF1* fitFunc = fitter.fit(hist, i, j, k, l);

        fillParams(fitFunc);

        downEdge = fitter.getDownEdge(hist, i, j, k, l);
        upEdge = fitter.getUpEdge(hist, i, j, k, l);
        width_mm = hist->GetBinWidth(1);
        
        // --- Расчет Yield и его погрешности ---
        // Используем std::abs, чтобы избежать проблем с отрицательными параметрами, если они есть
        yield = std::sqrt(TMath::Pi() / 2.0) * ampGauss * std::abs(stdDevGauss) * (1.0 + std::abs(asymGauss)) / width_mm;
        
        // Погрешность рассчитываем только если параметры не нулевые, чтобы избежать деления на ноль
        if (ampGauss != 0 && stdDevGauss != 0) {
            eYield = yield * std::sqrt(
                std::pow(eAmpGauss / ampGauss, 2) +
                std::pow(eStdDevGauss / stdDevGauss, 2) +
                std::pow(eAsymGauss / (1.0 + std::abs(asymGauss)), 2)
            );
        } else {
            eYield = 0.0;
        }
        // ---------------------------------------


        fitter.getNeutronIntegral(hist, i, j, k, l, neutronIntegral);
        fitter.getNeutronPeak(hist, i, j, k, l, neutronPosition, neutronValue);
        fitter.getMaxBin(hist, maxPosition, maxValue);

        index_q2 = i;
        index_w = j;
        index_cos_theta = k;
        index_phi = l;

        outputTree->Fill();
        delete fitFunc;
    }

    bool check(TH1F* hist) override {
        if (hist->GetEntries() > 0) return true;
        ++skippedHists;
        return false;
    }


    void logProgress() override {
        ++processedHists;
        if (processedHists % printEvery == 0 || processedHists == totalObjects) {
            log("Processing hist " + std::to_string(processedHists) + " / " + std::to_string(totalObjects), LogLevel::Info);
        }
    }

    bool finalize() override {
        outputFile->cd();
        outputTree->Write();

        log("Histograms fitting complete.", LogLevel::Info);
        log("Processed histogramms: " + std::to_string(processedHists), LogLevel::Info);
        log("Skipped histogramms due to empty or bad name: " + std::to_string(skippedHists), LogLevel::Info);

        return RootListObjectAnalysisStep<TH1F>::finalize();
    }

private:
    TTree* outputTree;
    HistogramFitterUniversal fitter;

    double ampGauss, meanGauss, stdDevGauss, asymGauss, ampGaussDelta, meanGaussDelta, stdDevGaussDelta, paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss, eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta, errA, errB, errC;
    double yield, eYield; // Новые переменные для выхода
    double downEdge, upEdge, width_mm, neutronIntegral, neutronValue, neutronPosition, maxValue, maxPosition;
    int index_q2, index_w, index_cos_theta, index_phi;

    void setupBranches() {
        outputTree->Branch("ampGauss", &ampGauss, "ampGauss/D");
        outputTree->Branch("meanGauss", &meanGauss, "meanGauss/D");
        outputTree->Branch("stdDevGauss", &stdDevGauss, "stdDevGauss/D");
        outputTree->Branch("asymGauss", &asymGauss, "asymGauss/D");
        outputTree->Branch("ampGaussDelta", &ampGaussDelta, "ampGaussDelta/D");
        outputTree->Branch("meanGaussDelta", &meanGaussDelta, "meanGaussDelta/D");
        outputTree->Branch("stdDevGaussDelta", &stdDevGaussDelta, "stdDevGaussDelta/D");
        outputTree->Branch("paramA", &paramA, "paramA/D");
        outputTree->Branch("paramB", &paramB, "paramB/D");
        outputTree->Branch("paramC", &paramC, "paramC/D");
        outputTree->Branch("eAmpGauss", &eAmpGauss, "eAmpGauss/D");
        outputTree->Branch("eMeanGauss", &eMeanGauss, "eMeanGauss/D");
        outputTree->Branch("eStdDevGauss", &eStdDevGauss, "eStdDevGauss/D");
        outputTree->Branch("eAsymGauss", &eAsymGauss, "eAsymGauss/D");
        outputTree->Branch("eAmpGaussDelta", &eAmpGaussDelta, "eAmpGaussDelta/D");
        outputTree->Branch("eMeanGaussDelta", &eMeanGaussDelta, "eMeanGaussDelta/D");
        outputTree->Branch("eStdDevGaussDelta", &eStdDevGaussDelta, "eStdDevGaussDelta/D");
        outputTree->Branch("errA", &errA, "errA/D");
        outputTree->Branch("errB", &errB, "errB/D");
        outputTree->Branch("errC", &errC, "errC/D");
        // Регистрация новых веток
        outputTree->Branch("yield", &yield, "yield/D");
        outputTree->Branch("eYield", &eYield, "eYield/D");

        outputTree->Branch("downEdge", &downEdge, "downEdge/D");
        outputTree->Branch("upEdge", &upEdge, "upEdge/D");
        outputTree->Branch("width_mm", &width_mm, "width_mm/D");
        outputTree->Branch("neutronIntegral", &neutronIntegral, "neutronIntegral/D");
        outputTree->Branch("neutronValue", &neutronValue, "neutronValue/D");
        outputTree->Branch("neutronPosition", &neutronPosition, "neutronPosition/D");
        outputTree->Branch("maxValue", &maxValue, "maxValue/D");

        outputTree->Branch("index_q2", &index_q2, "index_q2/I");
        outputTree->Branch("index_w", &index_w, "index_w/I");
        outputTree->Branch("index_cos_theta", &index_cos_theta, "index_cos_theta/I");
        outputTree->Branch("index_phi", &index_phi, "index_phi/I");
    }

    void fillParams(TF1* fitFunc) {
        ampGauss = fitFunc->GetParameter(0);
        meanGauss = fitFunc->GetParameter(1);
        stdDevGauss = fitFunc->GetParameter(2);
        asymGauss = fitFunc->GetParameter(3);
        ampGaussDelta = fitFunc->GetParameter(7);
        meanGaussDelta = fitFunc->GetParameter(8);
        stdDevGaussDelta = fitFunc->GetParameter(9);
        paramA = fitFunc->GetParameter(4);
        paramB = fitFunc->GetParameter(5);
        paramC = fitFunc->GetParameter(6);
        eAmpGauss = fitFunc->GetParError(0);
        eMeanGauss = fitFunc->GetParError(1);
        eStdDevGauss = fitFunc->GetParError(2);
        eAsymGauss = fitFunc->GetParError(3);
        eAmpGaussDelta = fitFunc->GetParError(7);
        eMeanGaussDelta = fitFunc->GetParError(8);
        eStdDevGaussDelta = fitFunc->GetParError(9);
        errA = fitFunc->GetParError(4);
        errB = fitFunc->GetParError(5);
        errC = fitFunc->GetParError(6);
    }
};
