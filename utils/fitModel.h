#pragma once

#include "TF1.h"
#include "TTree.h"
#include "TPaveText.h"

#include "legacyNeutronFitter.h"

class FitModel {
public:
    virtual ~FitModel() = default;

    // Фабричный метод для создания TF1 нужного типа
    virtual TF1* createTF1(const char* name, double xMin, double xMax) const = 0;
    
    // Функция фона
    virtual TF1* getBackgroundFormula(const char* name, double xMin, double xMax) const = 0;

    // Привязка веток TTree к внутренним буферам модели
    virtual void setupBranches(TTree* tree) = 0;

    // Установка параметров в функции перед отрисовкой
    virtual void setParams(TF1* fitFunc, TF1* bgFunc) const = 0;
};

class LegacyFitModel : public FitModel {
private:
    // Буферы для TTree (имена как в LegacyParamFitterStep)
    double amp, mean, sigma, asym;
    double pA, pB, pC;
    double ampD, meanD, sigmaD;

public:
    TF1* createTF1(const char* name, double xMin, double xMax) const override {
        std::string formula = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " 
                            "TMath::Max(0., [4]*(x - [5])*(x - [6])) + "
                            "[7]*exp(-0.5*((x-[8])/[9])^2)";
        return new TF1(name, formula.c_str(), xMin, xMax);
    }

    TF1* getBackgroundFormula(const char* name, double xMin, double xMax) const override {
        std::string formula = "TMath::Max(0., [0]*(x - [1])*(x - [2]))";
        return new TF1(name, formula.c_str(), xMin, xMax);
    }

    void setupBranches(TTree* tree) override {
        tree->SetBranchAddress("ampGauss", &amp);
        tree->SetBranchAddress("meanGauss", &mean);
        tree->SetBranchAddress("stdDevGauss", &sigma);
        tree->SetBranchAddress("asymGauss", &asym);

        tree->SetBranchAddress("paramA", &pA);
        tree->SetBranchAddress("paramB", &pB);
        tree->SetBranchAddress("paramC", &pC);
        
        tree->SetBranchAddress("ampGaussDelta", &ampD);
        tree->SetBranchAddress("meanGaussDelta", &meanD);
        tree->SetBranchAddress("stdDevGaussDelta", &sigmaD);
    }

    void setParams(TF1* f, TF1* b) const override {
        double pf[10] = {amp, mean, sigma, asym, pA, pB, pC, ampD, meanD, sigmaD};

        f->SetParameters(pf);
        if (b) {
            double pb[3] = {pA, pB, pC};
            b->SetParameters(pb);
        }
    }
};

class EMGFitModel : public FitModel {
private:
    // Буферы для веток TTree (специфичные для EMG)
    double yield, mu, sigma, tau;
    double ampBG, x1BG, x2BG;
    double ampD, meanD, sigmaD;

public:
    TF1* createTF1(const char* name, double xMin, double xMax) const override {
        return new TF1(name, EMGNeutronFitter::evaluateFullModel, xMin, xMax, 10);
    }

    TF1* getBackgroundFormula(const char* name, double xMin, double xMax) const override {
        // Фон у нас остался прежним - парабола
        std::string formula = "TMath::Max(0., [0]*(x - [1])*(x - [2]))";
        return new TF1(name, formula.c_str(), xMin, xMax);
    }

    void setupBranches(TTree* tree) override {
        // Здесь мы используем НОВЫЕ имена веток из EMGParamFitterStep
        tree->SetBranchAddress("yieldEMG", &yield);
        tree->SetBranchAddress("meanEMG", &mu);
        tree->SetBranchAddress("sigmaEMG", &sigma);
        tree->SetBranchAddress("tauEMG", &tau);
        
        tree->SetBranchAddress("ampBG", &ampBG);
        tree->SetBranchAddress("x1BG", &x1BG);
        tree->SetBranchAddress("x2BG", &x2BG);
        
        tree->SetBranchAddress("ampGaussDelta", &ampD);
        tree->SetBranchAddress("meanGaussDelta", &meanD);
        tree->SetBranchAddress("stdDevGaussDelta", &sigmaD);
    }

    void setParams(TF1* f, TF1* b) const override {
        double pf[10] = {yield, mu, sigma, tau, ampBG, x1BG, x2BG, ampD, meanD, sigmaD};
        f->SetParameters(pf);
        
        if (b) {
            double pb[3] = {ampBG, x1BG, x2BG};
            b->SetParameters(pb);
        }
    }
};