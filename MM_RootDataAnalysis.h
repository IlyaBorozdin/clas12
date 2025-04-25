#pragma once

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

#include "MM_project_const.h"
#include "n_piPlus_from_root.h"

using namespace std;

class MM_RootDataAnalysis {
public:
    MM_RootDataAnalysis() : file(MM_FILE.c_str(), "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << MM_FILE << std::endl;
            return;
        }

        missingMasses.resize(NUMBER_Q2, vector<TH1F*>(NUMBER_W, nullptr));
    }

    virtual ~MM_RootDataAnalysis() {
        if (file.IsOpen()) {
            file.Close();
        }
    }

    virtual void analysisEvent(int i, int j) = 0;
    virtual void results() = 0;

    void analysisCycle() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                missingMasses[i][j] = dynamic_cast<TH1F*>(file.Get(nameHistMM(i, j).c_str()));
                analysisEvent(i, j);
            }
        }
        results();
    }

protected:
    TFile file;
    vector<vector<TH1F*>> missingMasses;
};

class MM_ShortFitAnalysis : public MM_RootDataAnalysis {
public:
    MM_ShortFitAnalysis() {
        canvasesHistMM.resize(NUMBER_Q2, vector<TCanvas*>(NUMBER_FIGURE, nullptr));
        yieldVsW.resize(NUMBER_Q2, nullptr);
        fitFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        subFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));

        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_FIGURE; j++) {
                std::string nameCanvasMM = "canvas_mm_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                canvasesHistMM[i][j] = new TCanvas(nameCanvasMM.c_str(), "Histograms", 7200, 2700);                
            }
            for (int j = 0; j < NUMBER_W; j++) {
                std::string nameSubFunc = "SubFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);

                std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
                std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2) + " + strSubFunc;

                fitFuncs[i][j] = new TF1(nameFitFunc(i, j).c_str(), strFitFunc.c_str(), getDownEdge(i, j), getUpEdge(i, j));
                subFuncs[i][j] = new TF1(nameSubFunc.c_str(), strSubFunc.c_str(), getDownEdge(i, j), getUpEdge(i, j));
            }
            std::string nameCanvasY = "canvas_yield_" + std::to_string(i + 1);
            yieldVsW[i] = new TCanvas(nameCanvasY.c_str(), "Histogram", 1200, 900);
        }
    }

    ~MM_ShortFitAnalysis() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_FIGURE; j++) {
                if (canvasesHistMM[i][j]) delete canvasesHistMM[i][j];
            }
            for (int j = 0; j < NUMBER_W; j++) {
                if (fitFuncs[i][j]) delete fitFuncs[i][j];
                if (subFuncs[i][j]) delete subFuncs[i][j];
            }
            if (yieldVsW[i]) delete yieldVsW[i];    
        }
    }

    void analysisEvent(int i, int j) override {
        if (j % GRAPH_PER_FIGURE == 0) {
            canvasesHistMM[i][j / GRAPH_PER_FIGURE]->Divide(STR_PER_FIGURE, ROW_PER_FIGURE);
        }
        canvasesHistMM[i][j / GRAPH_PER_FIGURE]->cd((j % GRAPH_PER_FIGURE) + 1);
        
        int maxBin = missingMasses[i][j]->GetMaximumBin();
        double maxValue = missingMasses[i][j]->GetBinContent(maxBin);
        double maxPosition = missingMasses[i][j]->GetBinCenter(maxBin);

        fitFuncs[i][j]->SetParameter(0, maxValue);
        fitFuncs[i][j]->SetParameter(1, maxPosition);
        fitFuncs[i][j]->SetParameter(2, 0.03);
        fitFuncs[i][j]->SetParameter(4, 0.00);
        fitFuncs[i][j]->SetParameter(5, 0.00);
        fitFuncs[i][j]->SetParameter(6, 0.00);
        
        missingMasses[i][j]->Fit(nameFitFunc(i, j).c_str(), "R");

        subFuncs[i][j]->SetParameter(4, fitFuncs[i][j]->GetParameter(4));
        subFuncs[i][j]->SetParameter(5, fitFuncs[i][j]->GetParameter(5));
        subFuncs[i][j]->SetParameter(6, fitFuncs[i][j]->GetParameter(6));

        missingMasses[i][j]->Draw();
        fitFuncs[i][j]->Draw("same");
        subFuncs[i][j]->SetLineColor(kBlue);
        subFuncs[i][j]->Draw("same");

        if (j % GRAPH_PER_FIGURE == GRAPH_PER_FIGURE - 1) {
            canvasesHistMM[i][j / GRAPH_PER_FIGURE]->Update();
            std::string figNameCanvasMM = IMG_BUFF + "MM_elements_" + std::to_string(i + 1) + "-" + std::to_string(j / GRAPH_PER_FIGURE + 1) + ".png";
            canvasesHistMM[i][j / GRAPH_PER_FIGURE]->SaveAs(figNameCanvasMM.c_str());
        }
    }

    void results() override {
        for (int i = 0; i < NUMBER_Q2; i++) {
            std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TGraphErrors* graph = new TGraphErrors(NUMBER_W);
            graph->SetTitle(title.c_str());
            double yPoints[NUMBER_W];
            double eyPoints[NUMBER_W];

            for (int j = 0; j < NUMBER_W; j++) {
                double ampGauss = fitFuncs[i][j]->GetParameter(0);
                double stdDevGauss = abs(fitFuncs[i][j]->GetParameter(2));
                double eAmpGauss = fitFuncs[i][j]->GetParError(0);
                double eStdDevGauss = abs(fitFuncs[i][j]->GetParError(2));

                yPoints[j] = sqrt(2 * M_PI) * ampGauss * stdDevGauss / MM_BIN_WIDTH;
                eyPoints[j] = yPoints[j] * sqrt(pow(eAmpGauss / ampGauss, 2) + pow(eStdDevGauss /stdDevGauss, 2));

                graph->SetPoint(j, W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoints[j]);
                graph->SetPointError(j, 0.5 * (W_MAX - W_MIN) / NUMBER_W, eyPoints[j]);
            }

            graph->SetLineWidth(2);
            graph->SetLineStyle(1);

            yieldVsW[i]->Divide(1, 1);
            yieldVsW[i]->cd(1);
            graph->SetLineColor(kBlue);
            graph->Draw("APC");  
            yieldVsW[i]->Update();
            std::string fig_yield_name = IMG_BUFF + "MM_yield_" + std::to_string(i + 1) + ".png";
            yieldVsW[i]->SaveAs(fig_yield_name.c_str());

            delete graph;
        }
    }

private:
    const int GRAPH_PER_FIGURE = 9;
    const int STR_PER_FIGURE = 3;
    const int ROW_PER_FIGURE = GRAPH_PER_FIGURE / STR_PER_FIGURE;
    const int NUMBER_FIGURE = NUMBER_W / GRAPH_PER_FIGURE;

    vector<vector<TCanvas*>> canvasesHistMM;
    vector<TCanvas*> yieldVsW;

    vector<vector<TF1*>> fitFuncs;
    vector<vector<TF1*>> subFuncs;

    double getUpEdge(int i, int j) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 1.00, 1.00, 1.00, 1.02, 1.04, 1.07, 1.09, 1.10, 1.11, 1.12, 1.12, 1.10, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        return MM_UP_EDGES[(int) (j * n / NUMBER_W)];
    }

    double getDownEdge(int i, int j) {
        return 0.83;
    }

    std::string nameFitFunc(int i, int j) {
        return "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
    }
};

class MM_FullFitAnalysis : public MM_RootDataAnalysis {
public:
    MM_FullFitAnalysis() {
        canvasesHistMM.resize(NUMBER_Q2, vector<TCanvas*>(NUMBER_FIGURE, nullptr));
        yieldVsW.resize(NUMBER_Q2, nullptr);
        fitFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        subFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        fullFitFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));

        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_FIGURE; j++) {
                std::string nameCanvasMM = "canvas_mm_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                canvasesHistMM[i][j] = new TCanvas(nameCanvasMM.c_str(), "Histograms", 7200, 2700);                
            }
            for (int j = 0; j < NUMBER_W; j++) {
                std::string nameSubFunc = "SubFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);

                std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
                std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
                std::string fullFuncStr = strFitFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

                fitFuncs[i][j] = new TF1(nameFitFunc(i, j).c_str(), strFitFunc.c_str(), getDownEdge(i, j), getUpEdge(i, j));
                subFuncs[i][j] = new TF1(nameSubFunc.c_str(), strSubFunc.c_str(), getDownEdge(i, j), getUpEdge(i, j));
                fullFitFuncs[i][j] = new TF1(nameFullFitFunc(i, j).c_str(), fullFuncStr.c_str(), getDownEdge(i, j), getUpEdge(i, j));
            }
            std::string nameCanvasY = "canvas_yield_" + std::to_string(i + 1);
            yieldVsW[i] = new TCanvas(nameCanvasY.c_str(), "Histogram", 1200, 900);
        }
    }

    ~MM_FullFitAnalysis() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_FIGURE; j++) {
                if (canvasesHistMM[i][j]) delete canvasesHistMM[i][j];
            }
            for (int j = 0; j < NUMBER_W; j++) {
                if (fitFuncs[i][j]) delete fitFuncs[i][j];
                if (subFuncs[i][j]) delete subFuncs[i][j];
                if (fullFitFuncs[i][j]) delete fullFitFuncs[i][j];
            }
            if (yieldVsW[i]) delete yieldVsW[i];    
        }
    }

    void analysisEvent(int i, int j) override {
        if (j % GRAPH_PER_FIGURE == 0) {
            canvasesHistMM[i][j / GRAPH_PER_FIGURE]->Divide(STR_PER_FIGURE, ROW_PER_FIGURE);
        }
        canvasesHistMM[i][j / GRAPH_PER_FIGURE]->cd((j % GRAPH_PER_FIGURE) + 1);
        
        int maxBin = missingMasses[i][j]->GetMaximumBin();
        double maxValue = missingMasses[i][j]->GetBinContent(maxBin);
        double maxPosition = missingMasses[i][j]->GetBinCenter(maxBin);

        fitFuncs[i][j]->SetParameter(0, maxValue);
        fitFuncs[i][j]->SetParameter(1, maxPosition);
        fitFuncs[i][j]->SetParameter(2, 0.03);
        fitFuncs[i][j]->SetParameter(3, 1.00);
        fitFuncs[i][j]->SetParameter(4, 0.00);
        fitFuncs[i][j]->SetParameter(5, 0.00);
        fitFuncs[i][j]->SetParameter(6, 0.00);

        fitFuncs[i][j]->SetParLimits(0, 0.5 * maxValue, 1.5 * maxValue);
        fitFuncs[i][j]->SetParLimits(1, 0.85, 1.15);
        fitFuncs[i][j]->SetParLimits(2, 0.0085, 0.085);
        fitFuncs[i][j]->SetParLimits(3, 0.10, 10.0);
        
        missingMasses[i][j]->Fit(nameFitFunc(i, j).c_str(), "R");

        for (int k = 0; k < 7; k++) {
            fullFitFuncs[i][j]->SetParameter(k, fitFuncs[i][j]->GetParameter(k));
        }
        fullFitFuncs[i][j]->SetParameter(7, 0.00);
        fullFitFuncs[i][j]->SetParameter(8, getDeltaPeak(i, j));
        fullFitFuncs[i][j]->SetParameter(9, 0.04);

        fullFitFuncs[i][j]->SetParLimits(0, 0.7 * fitFuncs[i][j]->GetParameter(0), 1.3 * fitFuncs[i][j]->GetParameter(0));
        fullFitFuncs[i][j]->SetParLimits(1, 0.8 * fitFuncs[i][j]->GetParameter(1), 1.2 * fitFuncs[i][j]->GetParameter(1));
        fullFitFuncs[i][j]->SetParLimits(2, 0.8 * fitFuncs[i][j]->GetParameter(2), 1.2 * fitFuncs[i][j]->GetParameter(2));
        fullFitFuncs[i][j]->SetParLimits(3, 0.8 * fitFuncs[i][j]->GetParameter(3), 1.2 * fitFuncs[i][j]->GetParameter(3));
        fullFitFuncs[i][j]->SetParLimits(7, 0.00, 0.2 * fitFuncs[i][j]->GetParameter(0));
        fullFitFuncs[i][j]->SetParLimits(8, 0.95, 1.26);
        fullFitFuncs[i][j]->SetParLimits(9, 0.01, 0.1);                        

        missingMasses[i][j]->Fit(nameFullFitFunc(i, j).c_str(), "R");

        for (int k = 4; k < 7; k++) {
            subFuncs[i][j]->SetParameter(k, fullFitFuncs[i][j]->GetParameter(k));
        }

        missingMasses[i][j]->Draw();
        fullFitFuncs[i][j]->Draw("same");
        subFuncs[i][j]->SetLineColor(kBlue);
        subFuncs[i][j]->Draw("same");

        if (j % GRAPH_PER_FIGURE == GRAPH_PER_FIGURE - 1) {
            canvasesHistMM[i][j / GRAPH_PER_FIGURE]->Update();
            std::string figNameCanvasMM = IMG_BUFF + "MM_elements_" + std::to_string(i + 1) + "-" + std::to_string(j / GRAPH_PER_FIGURE + 1) + ".png";
            canvasesHistMM[i][j / GRAPH_PER_FIGURE]->SaveAs(figNameCanvasMM.c_str());
        }
    }

    void results() override {
        for (int i = 0; i < NUMBER_Q2; i++) {
            std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TGraphErrors* graph = new TGraphErrors(NUMBER_W);
            graph->SetTitle(title.c_str());
            double yPoints[NUMBER_W];
            double eyPoints[NUMBER_W];

            for (int j = 0; j < NUMBER_W; j++) {
                double ampGauss = fullFitFuncs[i][j]->GetParameter(0);
                double stdDevGauss = abs(fullFitFuncs[i][j]->GetParameter(2));
                double asymGauss = abs(fullFitFuncs[i][j]->GetParameter(3));
                double eAmpGauss = fullFitFuncs[i][j]->GetParError(0);
                double eStdDevGauss = fullFitFuncs[i][j]->GetParError(2);
                double eAsymGauss = fullFitFuncs[i][j]->GetParError(3);

                yPoints[j] = sqrt(M_PI / 2) * ampGauss * stdDevGauss * (1 + asymGauss) / MM_BIN_WIDTH;
                eyPoints[j] = yPoints[j] * sqrt(pow(eAmpGauss / ampGauss, 2) + pow(eStdDevGauss / stdDevGauss, 2) + pow(eAsymGauss / (1 + asymGauss), 2));

                graph->SetPoint(j, W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoints[j]);
                graph->SetPointError(j, 0.5 * (W_MAX - W_MIN) / NUMBER_W, eyPoints[j]);
            }

            graph->SetLineWidth(2);
            graph->SetLineStyle(1);

            yieldVsW[i]->Divide(1, 1);
            yieldVsW[i]->cd(1);
            graph->SetLineColor(kBlue);
            graph->Draw("APC");  
            yieldVsW[i]->Update();
            std::string fig_yield_name = IMG_BUFF + "MM_yield_" + std::to_string(i + 1) + ".png";
            yieldVsW[i]->SaveAs(fig_yield_name.c_str());

            delete graph;
        }
    }

private:
    const int GRAPH_PER_FIGURE = 9;
    const int STR_PER_FIGURE = 3;
    const int ROW_PER_FIGURE = GRAPH_PER_FIGURE / STR_PER_FIGURE;
    const int NUMBER_FIGURE = NUMBER_W / GRAPH_PER_FIGURE;

    vector<vector<TCanvas*>> canvasesHistMM;
    vector<TCanvas*> yieldVsW;

    vector<vector<TF1*>> fitFuncs;
    vector<vector<TF1*>> subFuncs;
    vector<vector<TF1*>> fullFitFuncs;

    double getUpEdge(int i, int j) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.14, 1.18, 1.23, 1.24, 1.24, 1.28, 1.30, 1.32, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34 };
        return MM_UP_EDGES[(int) (j * n / NUMBER_W)];
    }

    double getDownEdge(int i, int j) {
        return 0.80;
    }

    double getDeltaPeak(int i, int j) {
        const int n = 36;
        const double DELTA_PEAK[n] = { 0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.96, 0.96, 0.98, 1.08, 1.08, 1.12, 1.14, 1.16, 1.16, 1.18, 1.20, 1.20, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22 };
        return DELTA_PEAK[(int) (j * n / NUMBER_W)];
    }

    std::string nameFitFunc(int i, int j) {
        return "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
    }

    std::string nameFullFitFunc(int i, int j) {
        return "FullFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
    }
};

class MM_RootDataAnalysis_4D {
public:
    MM_RootDataAnalysis_4D() : file(MM_FILE.c_str(), "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << MM_FILE << std::endl;
            return;
        }

        std::cout << "RUN C1" << std::endl;
        missingMasses.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            missingMasses[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                missingMasses[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    missingMasses[i][j][k].resize(NUMBER_PHI, nullptr);
                }
                std::cout << "i = " << i << "; j = " << j << std::endl;
            }
        }
    }

    virtual ~MM_RootDataAnalysis_4D() {
        if (file.IsOpen()) {
            file.Close();
        }
    }

    virtual void analysisEvent(int i, int j, int k, int l) = 0;
    virtual void results() = 0;

    void analysisCycle() {
        std::cout << "RUN F" << std::endl;
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        missingMasses[i][j][k][l] = dynamic_cast<TH1F*>(file.Get(nameHist4D_MM(i, j, k, l).c_str()));
                        analysisEvent(i, j, k, l);
                    }
                }
                std::cout << "i = " << i << "; j = " << j << std::endl;
            }
        }
        results();
    }

protected:
    TFile file;
    vector<vector<vector<vector<TH1F*>>>> missingMasses;
};

class MM_EntriesAnalysis : public MM_RootDataAnalysis_4D {
public:
    MM_EntriesAnalysis() {
        h_entries_distribution = new TH1F("h_entries_distribution", "Histogramms' Number of Entries;Number of entries;Counts", 50, 0, 10000);
        cEntries = new TCanvas("c_entries", "Number of entries", 800, 600);
        
    }

    ~MM_EntriesAnalysis() {
        if (h_entries_distribution) delete h_entries_distribution;
        if (cEntries) delete cEntries;
    }

    void analysisEvent(int i, int j, int k, int l) override {
        h_entries_distribution->Fill(missingMasses[i][j][k][l]->GetEntries());
    }

    void results() override {
        cEntries->Divide(1, 1);
        cEntries->cd(1);
        h_entries_distribution->Draw();
        cEntries->Update();
        cEntries->SaveAs((IMG_BUFF + "Number_of_Entries.png").c_str());
    }

private:
    TH1F* h_entries_distribution;
    TCanvas* cEntries;
};

class MM_FullFitAnalysis_4D : public MM_RootDataAnalysis_4D {
public:
    MM_FullFitAnalysis_4D(): file("fitMM_25.root", "RECREATE"), tree("FitData", "Tree with peaks' parameters") {

        std::cout << "RUN C2" << std::endl;
        fitFuncs.resize(NUMBER_Q2);
        fullFitFuncs.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            fitFuncs[i].resize(NUMBER_W);
            fullFitFuncs[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                fitFuncs[i][j].resize(NUMBER_COS_THETA);
                fullFitFuncs[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    fitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    fullFitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
                        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
                        std::string fullFuncStr = strFitFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

                        fitFuncs[i][j][k][l] = new TF1(nameFitFunc(i, j, k, l).c_str(), strFitFunc.c_str(), getDownEdge(i, j, k, l), getUpEdge(i, j, k, l));
                        fullFitFuncs[i][j][k][l] = new TF1(nameFullFitFunc(i, j, k, l).c_str(), fullFuncStr.c_str(), getDownEdge(i, j, k, l), getUpEdge(i, j, k, l));
                    }
                }
                std::cout << "i = " << i << "; j = " << j << std::endl;
            }
        }

        tree.Branch("ampGauss", &ampGauss, "ampGauss/D");
        tree.Branch("meanGauss", &meanGauss, "meanGauss/D");
        tree.Branch("stdDevGauss", &stdDevGauss, "stdDevGauss/D");
        tree.Branch("asymGauss", &asymGauss, "asymGauss/D");
        tree.Branch("ampGaussDelta", &ampGaussDelta, "ampGaussDelta/D");
        tree.Branch("meanGaussDelta", &meanGaussDelta, "meanGaussDelta/D");
        tree.Branch("stdDevGaussDelta", &stdDevGaussDelta, "stdDevGaussDelta/D");
        tree.Branch("paramA", &paramA, "paramA/D");
        tree.Branch("paramB", &paramB, "paramB/D");
        tree.Branch("paramC", &paramC, "paramC/D");
        tree.Branch("eAmpGauss", &eAmpGauss, "eAmpGauss/D");
        tree.Branch("eMeanGauss", &eMeanGauss, "eMeanGauss/D");
        tree.Branch("eStdDevGauss", &eStdDevGauss, "eStdDevGauss/D");
        tree.Branch("eAsymGauss", &eAsymGauss, "eAsymGauss/D");
        tree.Branch("eAmpGaussDelta", &eAmpGaussDelta, "eAmpGaussDelta/D");
        tree.Branch("eMeanGaussDelta", &eMeanGaussDelta, "eMeanGaussDelta/D");
        tree.Branch("eStdDevGaussDelta", &eStdDevGaussDelta, "eStdDevGaussDelta/D");
        tree.Branch("errA", &errA, "errA/D");
        tree.Branch("errB", &errB, "errB/D");
        tree.Branch("errC", &errC, "errC/D");
        tree.Branch("index_w", &index_w, "index_w/I");
        tree.Branch("index_q2", &index_q2, "index_q2/I");
        tree.Branch("index_cos_theta", &index_cos_theta, "index_cos_theta/I");
        tree.Branch("index_phi", &index_phi, "index_phi/I");
    }

    ~MM_FullFitAnalysis_4D() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        if (fitFuncs[i][j][k][l]) delete fitFuncs[i][j][k][l];
                        if (fullFitFuncs[i][j][k][l]) delete fullFitFuncs[i][j][k][l];
                    }
                }
            }
        }
    }

    void analysisEvent(int i, int j, int k, int l) override {   
        if (missingMasses[i][j][k][l]->GetEntries() < ENTRIES_MIN) {
            return;
        }

        int bin_min = missingMasses[i][j][k][l]->FindBin(0.85);
        int bin_max = missingMasses[i][j][k][l]->FindBin(1.05);

        int maxBin = bin_min;
        double maxValue = missingMasses[i][j][k][l]->GetBinContent(bin_min);

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            double value = missingMasses[i][j][k][l]->GetBinContent(bin);
            if (value > maxValue) {
                maxValue = value;
                maxBin = bin;
            }
        }

        double maxPosition = missingMasses[i][j][k][l]->GetBinCenter(maxBin);

        fitFuncs[i][j][k][l]->SetParameter(0, maxValue);
        fitFuncs[i][j][k][l]->SetParameter(1, maxPosition);
        fitFuncs[i][j][k][l]->SetParameter(2, 0.03);
        fitFuncs[i][j][k][l]->SetParameter(3, 1.00);
        fitFuncs[i][j][k][l]->SetParameter(4, 0.00);
        fitFuncs[i][j][k][l]->SetParameter(5, 0.00);
        fitFuncs[i][j][k][l]->SetParameter(6, 0.00);

        fitFuncs[i][j][k][l]->SetParLimits(0, 0.1 * maxValue, maxValue);
        fitFuncs[i][j][k][l]->SetParLimits(1, 0.85, 1.05);
        fitFuncs[i][j][k][l]->SetParLimits(2, 0.0085, 0.085);
        fitFuncs[i][j][k][l]->SetParLimits(3, 0.10, 10.0);
        
        missingMasses[i][j][k][l]->Fit(nameFitFunc(i, j, k, l).c_str(), "R");

        for (int f = 0; f < 7; f++) {
            fullFitFuncs[i][j][k][l]->SetParameter(f, fitFuncs[i][j][k][l]->GetParameter(f));
        }
        fullFitFuncs[i][j][k][l]->SetParameter(7, 0.00);
        fullFitFuncs[i][j][k][l]->SetParameter(8, getDeltaPeak(i, j, k, l));
        fullFitFuncs[i][j][k][l]->SetParameter(9, 0.04);

        fullFitFuncs[i][j][k][l]->SetParLimits(0, 0.7 * fitFuncs[i][j][k][l]->GetParameter(0), 1.3 * fitFuncs[i][j][k][l]->GetParameter(0));
        fullFitFuncs[i][j][k][l]->SetParLimits(1, 0.8 * fitFuncs[i][j][k][l]->GetParameter(1), 1.2 * fitFuncs[i][j][k][l]->GetParameter(1));
        fullFitFuncs[i][j][k][l]->SetParLimits(2, 0.8 * fitFuncs[i][j][k][l]->GetParameter(2), 1.2 * fitFuncs[i][j][k][l]->GetParameter(2));
        fullFitFuncs[i][j][k][l]->SetParLimits(3, 0.8 * fitFuncs[i][j][k][l]->GetParameter(3), 1.2 * fitFuncs[i][j][k][l]->GetParameter(3));
        fullFitFuncs[i][j][k][l]->SetParLimits(7, 0.00, 1.2 * missingMasses[i][j][k][l]->GetBinContent(missingMasses[i][j][k][l]->FindBin(getDeltaPeak(i, j, k, l))));
        fullFitFuncs[i][j][k][l]->SetParLimits(8, 0.95, 1.26);
        fullFitFuncs[i][j][k][l]->SetParLimits(9, 0.01, 0.1);                        

        missingMasses[i][j][k][l]->Fit(nameFullFitFunc(i, j, k, l).c_str(), "R");

        ampGauss = fullFitFuncs[i][j][k][l]->GetParameter(0);
        meanGauss = fullFitFuncs[i][j][k][l]->GetParameter(1);
        stdDevGauss = fullFitFuncs[i][j][k][l]->GetParameter(2);
        asymGauss = fullFitFuncs[i][j][k][l]->GetParameter(3);
        ampGaussDelta = fullFitFuncs[i][j][k][l]->GetParameter(7);
        meanGaussDelta = fullFitFuncs[i][j][k][l]->GetParameter(8);
        stdDevGaussDelta = fullFitFuncs[i][j][k][l]->GetParameter(9);
        paramA = fullFitFuncs[i][j][k][l]->GetParameter(4);
        paramB = fullFitFuncs[i][j][k][l]->GetParameter(5);
        paramC = fullFitFuncs[i][j][k][l]->GetParameter(6);
        eAmpGauss = fullFitFuncs[i][j][k][l]->GetParError(0);
        eMeanGauss = fullFitFuncs[i][j][k][l]->GetParError(1);
        eStdDevGauss = fullFitFuncs[i][j][k][l]->GetParError(2);
        eAsymGauss = fullFitFuncs[i][j][k][l]->GetParError(3);
        eAmpGaussDelta = fullFitFuncs[i][j][k][l]->GetParError(7);
        eMeanGaussDelta = fullFitFuncs[i][j][k][l]->GetParError(8);
        eStdDevGaussDelta = fullFitFuncs[i][j][k][l]->GetParError(9);
        errA = fullFitFuncs[i][j][k][l]->GetParError(4);
        errB = fullFitFuncs[i][j][k][l]->GetParError(5);
        errC = fullFitFuncs[i][j][k][l]->GetParError(6);
        index_q2 = i;
        index_w = j;
        index_cos_theta = k;
        index_phi = l;

        tree.Fill();
    }

    void results() override {
        file.cd();
        tree.Write();
        file.Close();
    }

private:
    vector<vector<vector<vector<TF1*>>>> fitFuncs;
    vector<vector<vector<vector<TF1*>>>> fullFitFuncs;

    TFile file;
    TTree tree;

    double ampGauss, meanGauss, stdDevGauss, asymGauss, ampGaussDelta, meanGaussDelta, stdDevGaussDelta, paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss, eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta, errA, errB, errC;
    int index_w, index_q2, index_cos_theta, index_phi;

    double getUpEdge(int i, int j, int k, int l) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.14, 1.18, 1.23, 1.24, 1.24, 1.28, 1.30, 1.32, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34 };
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

    std::string nameFitFunc(int i, int j, int k, int l) {
        return "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }

    std::string nameFullFitFunc(int i, int j, int k, int l) {
        return "FullFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }
};