#pragma once

#include "paramFitterStepBase.h"
#include "utils/legacyNeutronFitter.h"
#include "utils/fitModel.h"

#include "utils/jsonParameterProvider.h"
#include "utils/baseHistogramAnalyzer.h"

class LegacyParamFitterStep : public ParamFitterStepBase {
private:
    // Специфичные переменные
    double ampGauss, meanGauss, stdDevGauss, asymGauss;
    double ampGaussDelta, meanGaussDelta, stdDevGaussDelta;
    double paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss;
    double eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta;
    double errA, errB, errC;

public:
    LegacyParamFitterStep(const std::string& stepName, 
                          const std::string& input, 
                          const std::string& output, 
                          const std::string& config)
        : ParamFitterStepBase(stepName, input, output, 
          std::make_unique<EMGNeutronFitter>(
            new JsonParameterProvider(config),
            new BaseHistogramAnalyzer()
          )) {}

    std::unique_ptr<FitModel> getFitModel() const {
        return std::make_unique<LegacyFitModel>();
    }

protected:
    void setupSpecificBranches() override {
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
    }

    void fillSpecificParams(TF1* fitFunc) override {
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

    void fillSpecificEmpty() override {
        ampGauss = 0.0;     eAmpGauss = 0.0;
        meanGauss = 0.9396;  eMeanGauss = 0.0;
        stdDevGauss = 0.03; eStdDevGauss = 0.0;
        asymGauss = 1.0;    eAsymGauss = 0.0;

        ampGaussDelta = 0.0; eAmpGaussDelta = 0.0;
        meanGaussDelta = 1.232; eMeanGaussDelta = 0.0; // Номинал Дельты
        stdDevGaussDelta = 0.03; eStdDevGaussDelta = 0.0;

        paramA = 0.0; errA = 0.0;
        paramB = 0.0; errB = 0.0;
        paramC = 0.0; errC = 0.0;
    }
};

class EMGParamFitterStep : public ParamFitterStepBase {
private:
    // Мы меняем названия переменных на более "физичные" для EMG
    double yieldEMG, meanEMG, sigmaEMG, tauEMG;
    double ampGaussDelta, meanGaussDelta, stdDevGaussDelta;
    double ampBG, x1BG, x2BG;
    
    // Ошибки
    double eYieldEMG, eMeanEMG, eSigmaEMG, eTauEMG;
    double eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta;
    double eAmpBG, eX1BG, eX2BG;

public:
    EMGParamFitterStep(const std::string& stepName, 
                       const std::string& input, 
                       const std::string& output, 
                       const std::string& config)
        : ParamFitterStepBase(stepName, input, output, 
          std::make_unique<EMGNeutronFitter>(
            new JsonParameterProvider(config),
            new BaseHistogramAnalyzer()
          )) {}

    std::unique_ptr<FitModel> getFitModel() const {
        return std::make_unique<EMGFitModel>();
    }

protected:
    void setupSpecificBranches() override {
        // Называем ветки так, чтобы в ROOT-файле было сразу понятно, что это EMG
        outputTree->Branch("yieldEMG", &yieldEMG, "yieldEMG/D");
        outputTree->Branch("meanEMG", &meanEMG, "meanEMG/D");
        outputTree->Branch("sigmaEMG", &sigmaEMG, "sigmaEMG/D");
        outputTree->Branch("tauEMG", &tauEMG, "tauEMG/D");
        
        outputTree->Branch("ampGaussDelta", &ampGaussDelta, "ampGaussDelta/D");
        outputTree->Branch("meanGaussDelta", &meanGaussDelta, "meanGaussDelta/D");
        outputTree->Branch("stdDevGaussDelta", &stdDevGaussDelta, "stdDevGaussDelta/D");
        
        outputTree->Branch("ampBG", &ampBG, "ampBG/D");
        outputTree->Branch("x1BG", &x1BG, "x1BG/D");
        outputTree->Branch("x2BG", &x2BG, "x2BG/D");

        // Ветки ошибок
        outputTree->Branch("eYieldEMG", &eYieldEMG, "eYieldEMG/D");
        outputTree->Branch("eMeanEMG", &eMeanEMG, "eMeanEMG/D");
        outputTree->Branch("eSigmaEMG", &eSigmaEMG, "eSigmaEMG/D");
        outputTree->Branch("eTauEMG", &eTauEMG, "eTauEMG/D");
        // ... и так далее для остальных ошибок
    }

    void fillSpecificParams(TF1* fitFunc) override {
        // Важно: в EMG параметр p[0] — это уже интеграл (Yield)!
        yieldEMG = fitFunc->GetParameter(0);
        meanEMG  = fitFunc->GetParameter(1);
        sigmaEMG = fitFunc->GetParameter(2);
        tauEMG   = fitFunc->GetParameter(3);
        
        ampBG = fitFunc->GetParameter(4);
        x1BG  = fitFunc->GetParameter(5);
        x2BG  = fitFunc->GetParameter(6);
        
        ampGaussDelta    = fitFunc->GetParameter(7);
        meanGaussDelta   = fitFunc->GetParameter(8);
        stdDevGaussDelta = fitFunc->GetParameter(9);

        // Ошибки
        eYieldEMG = fitFunc->GetParError(0);
        eMeanEMG  = fitFunc->GetParError(1);
        eSigmaEMG = fitFunc->GetParError(2);
        eTauEMG   = fitFunc->GetParError(3);

        eAmpBG = fitFunc->GetParError(4);
        eX1BG  = fitFunc->GetParError(5);
        eX2BG  = fitFunc->GetParError(6);

        eAmpGaussDelta    = fitFunc->GetParError(7);
        eMeanGaussDelta   = fitFunc->GetParError(8);
        eStdDevGaussDelta = fitFunc->GetParError(9);
    }

    void fillSpecificEmpty() override {
        yieldEMG = 0.0;     eYieldEMG = 0.0;
        meanEMG = 0.9396;    eMeanEMG = 0.0;
        sigmaEMG = 0.03;    eSigmaEMG = 0.0;
        tauEMG = 0.0;       eTauEMG = 0.0;

        ampGaussDelta = 0.0; 
        meanGaussDelta = 1.232; 
        stdDevGaussDelta = 0.03;

        ampBG = 0.0; x1BG = 0.0; x2BG = 0.0;
    }
};