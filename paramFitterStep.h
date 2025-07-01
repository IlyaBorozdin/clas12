#pragma once

#include "TH1F.h"
#include "TF1.h"
#include "TTree.h"

#include "MM_project_utils.h"
#include "source/ManageClasses/rootListObjectAnalysisStep.h"

class ParamFitterStep : public RootListObjectAnalysisStep<TH1F> {
public:
    ParamFitterStep(const std::string& stepName,
                    const std::string& inputFileName,
                    const std::string& outputFileName)
        : RootListObjectAnalysisStep<TH1F>(stepName, { {"main", inputFileName} }, cellGenerator, outputFileName) {}

    virtual ~ParamFitterStep() {
        delete outputTree;
    }

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
        if (!parseCellName(hist->GetName(), i, j, k, l)) {
            ++skippedHists;
            return;
        }

        TF1* fitFunc = histFit(hist, i, j, k, l);

        fillParams(fitFunc);
        downEdge = getDownEdge(i, j, k, l);
        upEdge = getUpEdge(i, j, k, l);
        width_mm = hist->GetBinWidth(1);
        index_q2 = i;
        index_w = j;
        index_cos_theta = k;
        index_phi = l;

        outputTree->Fill();
        delete fitFunc;
    }

    bool check(TH1F* hist) override {
        if (hist->GetEntries() > 0) {
            return true;
        }
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

        return RootListObjectAnalysisStep<TH1F>::finalize();
    }

private:
    TTree* outputTree;

    double ampGauss, meanGauss, stdDevGauss, asymGauss, ampGaussDelta, meanGaussDelta, stdDevGaussDelta, paramA, paramB, paramC, width_mm;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss, eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta, errA, errB, errC;
    double downEdge, upEdge;
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

        outputTree->Branch("downEdge", &downEdge, "downEdge/D");
        outputTree->Branch("upEdge", &upEdge, "upEdge/D");
        outputTree->Branch("width_mm", &width_mm, "width_mm/D");
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

    double getUpEdge(int i, int j, int k, int l) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.14, 1.18, 1.23, 1.24, 1.24, 1.28, 1.30, 1.32, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34 };
        return MM_UP_EDGES[(int) (j * n / NUMBER_W)];
    }

    double getUpShortEdge(int i, int j, int k, int l) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.09, 1.10, 1.11, 1.12, 1.12, 1.10, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        return MM_UP_EDGES[(int) (j * n / NUMBER_W)];
    }

    double getDownEdge(int i, int j, int k, int l) {
        return 0.80;
    }

    double getDeltaPeak(int i, int j, int k, int l) {
        const int n = 36;
        const double DELTA_PEAK[n] = { 0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.96, 0.96, 0.98, 1.08, 1.08, 1.12, 1.14, 1.16, 1.16, 1.18, 1.20, 1.20, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22 };
        return DELTA_PEAK[(int) (j * n / NUMBER_W)];
    }

    void getNeutronPeak(TH1F* hist, double& maxPosition, double& maxValue) {
        int bin_min = hist->FindBin(0.86);
        int bin_max = hist->FindBin(1.00);

        int maxBin = bin_min;
        maxValue = hist->GetBinContent(bin_min);

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            double value = hist->GetBinContent(bin);
            if (value > maxValue) {
                maxValue = value;
                maxBin = bin;
            }
        }

        maxPosition = hist->GetBinCenter(maxBin);
    }

    TF1* histFit(TH1F* hist, int i, int j, int k, int l) {
        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
        std::string strPeakFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
        std::string strFitFunc = strPeakFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        double upEdge = getUpEdge(i, j, k, l);
        double upShortEdge = getUpShortEdge(i, j, k, l);
        double downEdge = getDownEdge(i, j, k, l);
        double deltaPeak = getDeltaPeak(i, j, k, l);
        double maxPosition, maxValue;
        getNeutronPeak(hist, maxPosition, maxValue);

        TF1* fitFunc = new TF1("Fit Function", strFitFunc.c_str(), downEdge, upEdge);
        TF1* primaryFitFunc = new TF1("Primary Fit Function", strPeakFunc.c_str(), downEdge, upShortEdge);

        primaryFitFunc->SetParameter(0, maxValue);
        primaryFitFunc->SetParameter(1, maxPosition);
        primaryFitFunc->SetParameter(2, 0.03);
        primaryFitFunc->SetParameter(3, 1.00);
        primaryFitFunc->SetParameter(4, 0.00);
        primaryFitFunc->SetParameter(5, 0.80);
        primaryFitFunc->SetParameter(6, 1.60);

        primaryFitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        primaryFitFunc->SetParLimits(1, 0.86, 1.00);
        primaryFitFunc->SetParLimits(2, 0.0085, 0.085);
        primaryFitFunc->SetParLimits(3, 0.33, 3.00);
        primaryFitFunc->SetParLimits(5, 0.70, 0.95);
        primaryFitFunc->SetParLimits(6, 0.70, 4.00);

        hist->Fit(primaryFitFunc, "RQ");

        for (int f = 0; f < 7; f++) {
            fitFunc->SetParameter(f, primaryFitFunc->GetParameter(f));
        }
        fitFunc->SetParameter(7, 1.00);
        fitFunc->SetParameter(8, deltaPeak);
        fitFunc->SetParameter(9, 0.03);

        fitFunc->SetParLimits(0, 0.7 * primaryFitFunc->GetParameter(0), 1.3 * primaryFitFunc->GetParameter(0));
        fitFunc->SetParLimits(1, 0.9 * primaryFitFunc->GetParameter(1), 1.1 * primaryFitFunc->GetParameter(1));
        fitFunc->SetParLimits(2, 0.9 * primaryFitFunc->GetParameter(2), 1.1 * primaryFitFunc->GetParameter(2));
        fitFunc->SetParLimits(3, 0.9 * primaryFitFunc->GetParameter(3), 1.1 * primaryFitFunc->GetParameter(3));
        fitFunc->SetParLimits(7, 0.00, hist->GetBinContent(hist->FindBin(deltaPeak)));
        fitFunc->SetParLimits(8, 0.8 * deltaPeak, 1.2 * deltaPeak);
        fitFunc->SetParLimits(9, 0.01, 0.1);                        

        hist->Fit(fitFunc, "RQ");

        delete primaryFitFunc;
        return fitFunc;
    }
};
