#pragma once

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

// #include "MM_project_const.h"
#include "n_piPlus_from_root.h"

using namespace std;

/**
 * Базовый класс
 * работа с 2д гистограммами
*/
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

/**
 * По заготовленной гистограмме 2д:
 * производит аппроксимацию одним симметричным гауссом с фоном,
 * печетает графики
*/
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

/**
 * По заготовленной гистограмме 2д:
 * производит аппроксимацию двумя гауссами (ассиметрично) с фоном,
 * печетает графики
*/
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

/**
 * Базовый класс
 * работа с 4d гистограммами
*/
class MM_RootDataAnalysis_4D {
public:
    MM_RootDataAnalysis_4D(std::string readFileName = MM_FILE) : file(readFileName.c_str(), "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << readFileName << std::endl;
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
        std::cout << "FINISHED" << std::endl;
    }

    virtual void analysisEvent(int i, int j, int k, int l) = 0;
    virtual bool cutFunction(int i, int j, int k, int l) = 0;
    virtual void results() = 0;

    void analysisCycle() {
        std::cout << "RUN F" << std::endl;
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        missingMasses[i][j][k][l] = dynamic_cast<TH1F*>(file.Get(nameHist4D_MM(i, j, k, l).c_str()));
                        if (cutFunction(i, j, k, l)) {
                            analysisEvent(i, j, k, l);
                        }
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

/**
 * Побочный класс
 * по заготовленной гистограмме считает количество событий, строит график
*/
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

/**
 * Переделать,
 * реализует аппроксимацию в две фазы
*/
class MM_TwoGaussPeakFit {
public:
    MM_TwoGaussPeakFit(TF1* fullFitFunc) : fullFitFunc(fullFitFunc) {
        fitFunc = new TF1("func_A", strFitFunc.c_str());
        // fullFitFunc = new TF1("func_B", strFullFunc.c_str());
    };

    ~MM_TwoGaussPeakFit() {
        if (fitFunc) delete fitFunc;
        // if (fullFitFunc) delete fullFitFunc;
    }

    void setShortEdges(double downEdge, double upEdge) {
        fitFunc->SetRange(downEdge, upEdge);
    }

    void setFullEdges(double downEdge, double upEdge) {
        fullFitFunc->SetRange(downEdge, upEdge);
    }

    void setDeltaPeak(double dp) {
        deltaPeak = dp;
    }

    // std::string getStrSubFunc() {
    //     return strSubFunc;
    // }

    std::string getStrFitFunc() {
        return strFitFunc;
    }

    // std::string getStrFullFunc() {
    //     return strFullFunc;
    // }

    void fitting(TH1F* missingMass) {
        int bin_min = missingMass->FindBin(0.86);
        int bin_max = missingMass->FindBin(1.00);

        int maxBin = bin_min;
        double maxValue = missingMass->GetBinContent(bin_min);

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            double value = missingMass->GetBinContent(bin);
            if (value > maxValue) {
                maxValue = value;
                maxBin = bin;
            }
        }

        double maxPosition = missingMass->GetBinCenter(maxBin);

        fitFunc->SetParameter(0, maxValue);
        fitFunc->SetParameter(1, maxPosition);
        fitFunc->SetParameter(2, 0.03);
        fitFunc->SetParameter(3, 1.00);
        fitFunc->SetParameter(4, 0.00);
        fitFunc->SetParameter(5, 0.00);
        fitFunc->SetParameter(6, 0.00);

        fitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        fitFunc->SetParLimits(1, 0.86, 1.00);
        fitFunc->SetParLimits(2, 0.0085, 0.085);
        fitFunc->SetParLimits(3, 0.33, 3.00);
        
        missingMass->Fit(fitFunc, "RQ");

        for (int f = 0; f < 7; f++) {
            fullFitFunc->SetParameter(f, fitFunc->GetParameter(f));
        }
        fullFitFunc->SetParameter(7, 1.00);
        fullFitFunc->SetParameter(8, deltaPeak);
        fullFitFunc->SetParameter(9, 0.03);

        fullFitFunc->SetParLimits(0, 0.7 * fitFunc->GetParameter(0), 1.3 * fitFunc->GetParameter(0));
        fullFitFunc->SetParLimits(1, 0.9 * fitFunc->GetParameter(1), 1.1 * fitFunc->GetParameter(1));
        fullFitFunc->SetParLimits(2, 0.9 * fitFunc->GetParameter(2), 1.1 * fitFunc->GetParameter(2));
        fullFitFunc->SetParLimits(3, 0.9 * fitFunc->GetParameter(3), 1.1 * fitFunc->GetParameter(3));
        fullFitFunc->SetParLimits(7, 0.00, missingMass->GetBinContent(missingMass->FindBin(deltaPeak)));
        fullFitFunc->SetParLimits(8, 0.8 * deltaPeak, 1.2 * deltaPeak);
        fullFitFunc->SetParLimits(9, 0.01, 0.1);                        

        missingMass->Fit(fullFitFunc, "RQ");
    }
private:
    // std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
    std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + [4]*x*x + [5]*x + [6]";
    // std::string strFullFunc = strFitFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

    TF1* fitFunc;
    TF1* fullFitFunc;

    double deltaPeak;
};

/**
 * По гистограммам 4д 
 * проводит аппроксимацию с использованием two_gauss_peak_fit
 * создаёт файл дерево с параметрами аппроксимации
*/
class MM_FullFitAnalysis_4D : public MM_RootDataAnalysis_4D {
public:
    MM_FullFitAnalysis_4D(const char* outFileName): file(outFileName, "RECREATE"), tree("FitData", "Tree with peaks' parameters") {
        std::cout << "RUN C2" << std::endl;
        // fitFuncs.resize(NUMBER_Q2);
        fullFitFuncs.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            // fitFuncs[i].resize(NUMBER_W);
            fullFitFuncs[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                // fitFuncs[i][j].resize(NUMBER_COS_THETA);
                fullFitFuncs[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    // fitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    fullFitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
                        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
                        std::string fullFuncStr = strFitFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

                        // fitFuncs[i][j][k][l] = new TF1(nameFitFunc(i, j, k, l).c_str(), strFitFunc.c_str(), getDownEdge(i, j, k, l), getUpShortEdge(i, j, k, l));
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
                        // if (fitFuncs[i][j][k][l]) delete fitFuncs[i][j][k][l];
                        if (fullFitFuncs[i][j][k][l]) delete fullFitFuncs[i][j][k][l];
                    }
                }
            }
        }
    }

    bool cutFunction(int i, int j, int k, int l) override {
        return !(!missingMasses[i][j][k][l] || missingMasses[i][j][k][l]->GetEntries() < ENTRIES_MIN);
    }

    void analysisEvent(int i, int j, int k, int l) override {   
        /*
        int bin_min = missingMasses[i][j][k][l]->FindBin(0.86);
        int bin_max = missingMasses[i][j][k][l]->FindBin(1.00);

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

        fitFuncs[i][j][k][l]->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        fitFuncs[i][j][k][l]->SetParLimits(1, 0.86, 1.00);
        fitFuncs[i][j][k][l]->SetParLimits(2, 0.0085, 0.085);
        fitFuncs[i][j][k][l]->SetParLimits(3, 0.33, 3.00);
        
        auto fitResult = missingMasses[i][j][k][l]->Fit(fitFuncs[i][j][k][l], "RQ");
        if (fitResult.Get() == nullptr) {
            std::cerr << "WARNING: First fit failed for i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << std::endl;
        }

        for (int f = 0; f < 7; f++) {
            fullFitFuncs[i][j][k][l]->SetParameter(f, fitFuncs[i][j][k][l]->GetParameter(f));
        }
        fullFitFuncs[i][j][k][l]->SetParameter(7, 1.00);
        fullFitFuncs[i][j][k][l]->SetParameter(8, getDeltaPeak(i, j, k, l));
        fullFitFuncs[i][j][k][l]->SetParameter(9, 0.03);

        fullFitFuncs[i][j][k][l]->SetParLimits(0, 0.7 * fitFuncs[i][j][k][l]->GetParameter(0), 1.3 * fitFuncs[i][j][k][l]->GetParameter(0));
        fullFitFuncs[i][j][k][l]->SetParLimits(1, 0.9 * fitFuncs[i][j][k][l]->GetParameter(1), 1.1 * fitFuncs[i][j][k][l]->GetParameter(1));
        fullFitFuncs[i][j][k][l]->SetParLimits(2, 0.9 * fitFuncs[i][j][k][l]->GetParameter(2), 1.1 * fitFuncs[i][j][k][l]->GetParameter(2));
        fullFitFuncs[i][j][k][l]->SetParLimits(3, 0.9 * fitFuncs[i][j][k][l]->GetParameter(3), 1.1 * fitFuncs[i][j][k][l]->GetParameter(3));
        fullFitFuncs[i][j][k][l]->SetParLimits(7, 0.00, missingMasses[i][j][k][l]->GetBinContent(missingMasses[i][j][k][l]->FindBin(getDeltaPeak(i, j, k, l))));
        fullFitFuncs[i][j][k][l]->SetParLimits(8, 0.8 * getDeltaPeak(i, j, k, l), 1.2 * getDeltaPeak(i, j, k, l));
        fullFitFuncs[i][j][k][l]->SetParLimits(9, 0.01, 0.1);                        

        fitResult = missingMasses[i][j][k][l]->Fit(fullFitFuncs[i][j][k][l], "RQ");
        if (fitResult.Get() == nullptr) {
            std::cerr << "WARNING: Second fit failed for i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << std::endl;
        }
        */
        MM_TwoGaussPeakFit fittingProcedure(fullFitFuncs[i][j][k][l]);
        fittingProcedure.setShortEdges(getDownEdge(i, j, k, l), getUpShortEdge(i, j, k, l));
        fittingProcedure.setDeltaPeak(getDeltaPeak(i, j, k, l));
        fittingProcedure.fitting(missingMasses[i][j][k][l]);

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
    // vector<vector<vector<vector<TF1*>>>> fitFuncs;
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

    // std::string nameFitFunc(int i, int j, int k, int l) {
    //     return "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    // }

    std::string nameFullFitFunc(int i, int j, int k, int l) {
        return "FullFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }
};

/**
 * Базовый класс
 * для работы с деревом параметров аппроксимации
 * рисует графики
*/
class MM_FittedDataAnalysis_4D {
public:
    MM_FittedDataAnalysis_4D(const char* outFileName) 
        : file(outFileName, "READ"), tree(nullptr), file_hists(MM_FILE.c_str(), "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << outFileName << std::endl;
            return;
        }

        if (!file_hists.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << MM_FILE << std::endl;
            return;
        }

        file.GetObject("FitData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'FitData' not found in file " << outFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("ampGauss", &ampGauss);
        tree->SetBranchAddress("meanGauss", &meanGauss);
        tree->SetBranchAddress("stdDevGauss", &stdDevGauss);
        tree->SetBranchAddress("asymGauss", &asymGauss);
        tree->SetBranchAddress("ampGaussDelta", &ampGaussDelta);
        tree->SetBranchAddress("meanGaussDelta", &meanGaussDelta);
        tree->SetBranchAddress("stdDevGaussDelta", &stdDevGaussDelta);
        tree->SetBranchAddress("paramA", &paramA);
        tree->SetBranchAddress("paramB", &paramB);
        tree->SetBranchAddress("paramC", &paramC);
        tree->SetBranchAddress("eAmpGauss", &eAmpGauss);
        tree->SetBranchAddress("eMeanGauss", &eMeanGauss);
        tree->SetBranchAddress("eStdDevGauss", &eStdDevGauss);
        tree->SetBranchAddress("eAsymGauss", &eAsymGauss);
        tree->SetBranchAddress("eAmpGaussDelta", &eAmpGaussDelta);
        tree->SetBranchAddress("eMeanGaussDelta", &eMeanGaussDelta);
        tree->SetBranchAddress("eStdDevGaussDelta", &eStdDevGaussDelta);
        tree->SetBranchAddress("errA", &errA);
        tree->SetBranchAddress("errB", &errB);
        tree->SetBranchAddress("errC", &errC);
        tree->SetBranchAddress("index_w", &index_w);
        tree->SetBranchAddress("index_q2", &index_q2);
        tree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        tree->SetBranchAddress("index_phi", &index_phi);
    }

    ~MM_FittedDataAnalysis_4D() {
        if (file.IsOpen()) {
            file.Close();
        }
        if (file_hists.IsOpen()) {
            file_hists.Close();
        }
    }

    bool cutFunction() {
        if (meanGauss < 0.86 || meanGauss > 1.00) return false;
        if (stdDevGauss < 0.0085 || stdDevGauss > 0.085) return false;
        if (asymGauss < 0.25 || asymGauss > 4.00) return false;
        if (ampGaussDelta < 0.0) return false;
        // if (meanGaussDelta < 0.956 || meanGaussDelta > 1.254) return false;
        // if (stdDevGaussDelta < 0.01 || stdDevGaussDelta > 0.10) return false;
        return true;
    }

    void analysisEvent() {
        TH1F* missingMass = dynamic_cast<TH1F*>(file_hists.Get(nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi).c_str()));
        TCanvas* canvas = new TCanvas("canvas", "Histogram", 1200, 900);

        std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        TF1* fitFuncs = new TF1("FitFunc", strFitFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));
        TF1* subFuncs = new TF1("SubFunc", strSubFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));

        fitFuncs->SetParameter(0, ampGauss);
        fitFuncs->SetParameter(1, meanGauss);
        fitFuncs->SetParameter(2, stdDevGauss);
        fitFuncs->SetParameter(3, asymGauss);
        fitFuncs->SetParameter(4, paramA);
        fitFuncs->SetParameter(5, paramB);
        fitFuncs->SetParameter(6, paramC);
        fitFuncs->SetParameter(7, ampGaussDelta);
        fitFuncs->SetParameter(8, meanGaussDelta);
        fitFuncs->SetParameter(9, stdDevGaussDelta);

        subFuncs->SetParameter(4, paramA);
        subFuncs->SetParameter(5, paramB);
        subFuncs->SetParameter(6, paramC);

        canvas->Divide(1, 1);
        canvas->cd(1);

        missingMass->Draw();
        fitFuncs->Draw("same");
        subFuncs->SetLineColor(kBlue);
        subFuncs->Draw("same");

        canvas->Update();
        std::string figNameCanvasMM = IMG_BUFF + "big50/" + nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi) + ".png";
        canvas->SaveAs(figNameCanvasMM.c_str());

        delete canvas;
        delete fitFuncs;
        delete subFuncs;
    }

    void results() {}

    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            if (cutFunction()) {
                analysisEvent();
            }
        }
        results();
    }

protected:
    TFile file;
    TTree* tree;

    TFile file_hists;

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
};

/**
 * Использовался однозначно для презентации
 * Просто подготавливает срез экспериментальных данных после фитирования
 * пробегается по делеву с параметрами аппроксимации
 * строит графики
*/
class MM_FittedDataAnalysis_4D_PRESENTATION {
public:
    MM_FittedDataAnalysis_4D_PRESENTATION(const char* outFileName) 
        : file(outFileName, "READ"), tree(nullptr), file_hists(MM_FILE.c_str(), "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << outFileName << std::endl;
            return;
        }

        if (!file_hists.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << MM_FILE << std::endl;
            return;
        }

        file.GetObject("FitData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'FitData' not found in file " << outFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("ampGauss", &ampGauss);
        tree->SetBranchAddress("meanGauss", &meanGauss);
        tree->SetBranchAddress("stdDevGauss", &stdDevGauss);
        tree->SetBranchAddress("asymGauss", &asymGauss);
        tree->SetBranchAddress("ampGaussDelta", &ampGaussDelta);
        tree->SetBranchAddress("meanGaussDelta", &meanGaussDelta);
        tree->SetBranchAddress("stdDevGaussDelta", &stdDevGaussDelta);
        tree->SetBranchAddress("paramA", &paramA);
        tree->SetBranchAddress("paramB", &paramB);
        tree->SetBranchAddress("paramC", &paramC);
        tree->SetBranchAddress("eAmpGauss", &eAmpGauss);
        tree->SetBranchAddress("eMeanGauss", &eMeanGauss);
        tree->SetBranchAddress("eStdDevGauss", &eStdDevGauss);
        tree->SetBranchAddress("eAsymGauss", &eAsymGauss);
        tree->SetBranchAddress("eAmpGaussDelta", &eAmpGaussDelta);
        tree->SetBranchAddress("eMeanGaussDelta", &eMeanGaussDelta);
        tree->SetBranchAddress("eStdDevGaussDelta", &eStdDevGaussDelta);
        tree->SetBranchAddress("errA", &errA);
        tree->SetBranchAddress("errB", &errB);
        tree->SetBranchAddress("errC", &errC);
        tree->SetBranchAddress("index_w", &index_w);
        tree->SetBranchAddress("index_q2", &index_q2);
        tree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        tree->SetBranchAddress("index_phi", &index_phi);

        canvasHist.resize(100, nullptr);
        canvasYeild.resize(100, nullptr);
        graphs.resize(4);
        for (auto& vec : graphs) {
            vec.resize(100);
            for (auto& g : vec) {
                g = new TGraphErrors(NUMBER_W);
            }
        }
        for (int i = 0; i < 100; i++) {
            canvasHist[i] = new TCanvas(("canvasH" + std::to_string(i + 1)).c_str(), "HistogramH", 4800, 3600);
            canvasHist[i]->Divide(4, 3);
            canvasYeild[i] = new TCanvas(("canvasY" + std::to_string(i + 1)).c_str(), "HistogramY", 1200, 900);
            canvasYeild[i]->Divide(2, 2);
            for (int j = 0; j < 4; j++) {
                std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * j], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * j + 1], 1) + " GeV^{2}, cell-" + std::to_string(i + 1) + "; W, GeV; Yield";
                graphs[j][i]->SetTitle(title.c_str());
            }
        }
    }

    ~MM_FittedDataAnalysis_4D_PRESENTATION() {
        if (file.IsOpen()) {
            file.Close();
        }
        if (file_hists.IsOpen()) {
            file_hists.Close();
        }

        for (int i = 0; i < 100; i++) {
            if (canvasHist[i]) delete canvasHist[i];
            if (canvasYeild[i]) delete canvasYeild[i];
            for (int j = 0; j < 4; j++) {
                if (graphs[j][i]) delete graphs[j][i];
            }
        }
    }

    bool cutFunction() {
        if (meanGauss < 0.86 || meanGauss > 1.00) return false;
        if (stdDevGauss < 0.0085 || stdDevGauss > 0.085) return false;
        if (asymGauss < 0.25 || asymGauss > 4.00) return false;
        if (ampGaussDelta < 0.0) return false;
        // if (meanGaussDelta < 0.956 || meanGaussDelta > 1.254) return false;
        // if (stdDevGaussDelta < 0.01 || stdDevGaussDelta > 0.10) return false;
        return true;
    }

    void analysisEvent() {
        TH1F* missingMass = dynamic_cast<TH1F*>(file_hists.Get(nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi).c_str()));

        std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        TF1* fitFuncs = new TF1("FitFunc", strFitFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));
        TF1* subFuncs = new TF1("SubFunc", strSubFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));

        fitFuncs->SetParameter(0, ampGauss);
        fitFuncs->SetParameter(1, meanGauss);
        fitFuncs->SetParameter(2, stdDevGauss);
        fitFuncs->SetParameter(3, asymGauss);
        fitFuncs->SetParameter(4, paramA);
        fitFuncs->SetParameter(5, paramB);
        fitFuncs->SetParameter(6, paramC);
        fitFuncs->SetParameter(7, ampGaussDelta);
        fitFuncs->SetParameter(8, meanGaussDelta);
        fitFuncs->SetParameter(9, stdDevGaussDelta);

        subFuncs->SetParameter(4, paramA);
        subFuncs->SetParameter(5, paramB);
        subFuncs->SetParameter(6, paramC);

        int numberRes = (index_w == 2) ? 1 : (index_w == 7) ? 2 : (index_w == 11) ? 3 : -1;
        if (numberRes > 0) {
            canvasHist[10 * index_cos_theta + index_phi]->cd(4 * (numberRes - 1) + index_q2 + 1);
            missingMass->Draw();
            fitFuncs->Draw("same");
            subFuncs->SetLineColor(kBlue);
            subFuncs->Draw("same");
        }
        
        double yPoint = sqrt(M_PI / 2) * ampGauss * abs(stdDevGauss) * (1 + abs(asymGauss)) / MM_BIN_WIDTH;
        double yErr = yPoint * sqrt(pow(eAmpGauss / ampGauss, 2) + pow(eStdDevGauss / stdDevGauss, 2) + pow(eAsymGauss / (1 + asymGauss), 2));

        graphs[index_q2][10 * index_cos_theta + index_phi]->SetPoint(index_w, W_MIN + (index_w + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoint);
        graphs[index_q2][10 * index_cos_theta + index_phi]->SetPointError(index_w, 0.5 * (W_MAX - W_MIN) / NUMBER_W, yErr);
    }

    void results() {
        for (int i = 0; i < 100; i++) {
            canvasHist[i]->Update();
            std::string figNameCanvasRes = IMG_BUFF + "presentation50/resonances_" + std::to_string(i + 1) +  + ".png";
            canvasHist[i]->SaveAs(figNameCanvasRes.c_str());

            for (int j = 0; j < 4; j++) {
                graphs[j][i]->SetLineWidth(2);
                graphs[j][i]->SetLineStyle(1);
                canvasYeild[i]->cd(j + 1);
                graphs[j][i]->SetLineColor(kBlue);
                graphs[j][i]->Draw("APC");
            }
            canvasYeild[i]->Update();
            std::string figNameCanvasYeild = IMG_BUFF + "presentation50/yeild_" + std::to_string(i + 1) +  + ".png";
            canvasYeild[i]->SaveAs(figNameCanvasYeild.c_str());
        }
    }

    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            if (cutFunction()) {
                analysisEvent();
            }
        }
        results();
    }

protected:
    TFile file;
    TTree* tree;

    TFile file_hists;

    vector<TCanvas*> canvasHist;
    vector<TCanvas*> canvasYeild;
    vector<vector<TGraphErrors*>> graphs;

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
};


/**
class MM_ShortFitAnalysis_4D : public MM_RootDataAnalysis_4D {
public:
    MM_ShortFitAnalysis_4D(const char* outFileName): file(outFileName, "RECREATE"), tree("FitData", "Tree with peaks' parameters") {
        std::cout << "RUN C2" << std::endl;
        fitFuncs.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            fitFuncs[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                fitFuncs[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    fitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2) + [4]*x*x + [5]*x + [6]";
                        fitFuncs[i][j][k][l] = new TF1(nameFitFunc(i, j, k, l).c_str(), strFitFunc.c_str(), getDownEdge(i, j, k, l), getUpEdge(i, j, k, l));
                    }
                }
                std::cout << "i = " << i << "; j = " << j << std::endl;
            }
        }

        tree.Branch("ampGauss", &ampGauss, "ampGauss/D");
        tree.Branch("meanGauss", &meanGauss, "meanGauss/D");
        tree.Branch("stdDevGauss", &stdDevGauss, "stdDevGauss/D");
        tree.Branch("paramA", &paramA, "paramA/D");
        tree.Branch("paramB", &paramB, "paramB/D");
        tree.Branch("paramC", &paramC, "paramC/D");
        tree.Branch("eAmpGauss", &eAmpGauss, "eAmpGauss/D");
        tree.Branch("eMeanGauss", &eMeanGauss, "eMeanGauss/D");
        tree.Branch("eStdDevGauss", &eStdDevGauss, "eStdDevGauss/D");
        tree.Branch("errA", &errA, "errA/D");
        tree.Branch("errB", &errB, "errB/D");
        tree.Branch("errC", &errC, "errC/D");
        tree.Branch("index_w", &index_w, "index_w/I");
        tree.Branch("index_q2", &index_q2, "index_q2/I");
        tree.Branch("index_cos_theta", &index_cos_theta, "index_cos_theta/I");
        tree.Branch("index_phi", &index_phi, "index_phi/I");
    }

    ~MM_ShortFitAnalysis_4D() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        if (fitFuncs[i][j][k][l]) delete fitFuncs[i][j][k][l];
                    }
                }
            }
        }
    }

    bool cutFunction(int i, int j, int k, int l) override {
        return !(!missingMasses[i][j][k][l] || missingMasses[i][j][k][l]->GetEntries() < ENTRIES_MIN);
    }

    void analysisEvent(int i, int j, int k, int l) override {   
        int bin_min = missingMasses[i][j][k][l]->FindBin(0.85);
        int bin_max = missingMasses[i][j][k][l]->FindBin(0.98);

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
        fitFuncs[i][j][k][l]->SetParameter(4, 0.00);
        fitFuncs[i][j][k][l]->SetParameter(5, 0.00);
        fitFuncs[i][j][k][l]->SetParameter(6, 0.00);

        fitFuncs[i][j][k][l]->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        fitFuncs[i][j][k][l]->SetParLimits(1, 0.80, 1.00);
        fitFuncs[i][j][k][l]->SetParLimits(2, 0.0085, 0.085);
        
        auto fitResult = missingMasses[i][j][k][l]->Fit(fitFuncs[i][j][k][l], "RQ");
        if (fitResult.Get() == nullptr) {
            std::cerr << "WARNING: Fit failed for i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << std::endl;
        }

        ampGauss = fitFuncs[i][j][k][l]->GetParameter(0);
        meanGauss = fitFuncs[i][j][k][l]->GetParameter(1);
        stdDevGauss = fitFuncs[i][j][k][l]->GetParameter(2);
        paramA = fitFuncs[i][j][k][l]->GetParameter(4);
        paramB = fitFuncs[i][j][k][l]->GetParameter(5);
        paramC = fitFuncs[i][j][k][l]->GetParameter(6);
        eAmpGauss = fitFuncs[i][j][k][l]->GetParError(0);
        eMeanGauss = fitFuncs[i][j][k][l]->GetParError(1);
        eStdDevGauss = fitFuncs[i][j][k][l]->GetParError(2);
        errA = fitFuncs[i][j][k][l]->GetParError(4);
        errB = fitFuncs[i][j][k][l]->GetParError(5);
        errC = fitFuncs[i][j][k][l]->GetParError(6);
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

    TFile file;
    TTree tree;

    double ampGauss, meanGauss, stdDevGauss, paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, errA, errB, errC;
    int index_w, index_q2, index_cos_theta, index_phi;

    double getUpEdge(int i, int j, int k, int l) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 1.00, 1.00, 1.00, 1.02, 1.04, 1.07, 1.09, 1.10, 1.11, 1.12, 1.12, 1.10, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        return MM_UP_EDGES[(int) (j * n / NUMBER_W)];
    }

    double getDownEdge(int i, int j, int k, int l) {
        return 0.83;
    }

    std::string nameFitFunc(int i, int j, int k, int l) {
        return "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }

    std::string nameFullFitFunc(int i, int j, int k, int l) {
        return "FullFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }
};

class MM_FittedShortDataAnalysis_4D {
public:
    MM_FittedShortDataAnalysis_4D(const char* outFileName) 
        : file(outFileName, "READ"), tree(nullptr), file_hists(MM_FILE.c_str(), "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << outFileName << std::endl;
            return;
        }

        if (!file_hists.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << MM_FILE << std::endl;
            return;
        }

        file.GetObject("FitData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'FitData' not found in file " << outFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("ampGauss", &ampGauss);
        tree->SetBranchAddress("meanGauss", &meanGauss);
        tree->SetBranchAddress("stdDevGauss", &stdDevGauss);
        tree->SetBranchAddress("paramA", &paramA);
        tree->SetBranchAddress("paramB", &paramB);
        tree->SetBranchAddress("paramC", &paramC);
        tree->SetBranchAddress("eAmpGauss", &eAmpGauss);
        tree->SetBranchAddress("eMeanGauss", &eMeanGauss);
        tree->SetBranchAddress("eStdDevGauss", &eStdDevGauss);
        tree->SetBranchAddress("errA", &errA);
        tree->SetBranchAddress("errB", &errB);
        tree->SetBranchAddress("errC", &errC);
        tree->SetBranchAddress("index_w", &index_w);
        tree->SetBranchAddress("index_q2", &index_q2);
        tree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        tree->SetBranchAddress("index_phi", &index_phi);
    }

    ~MM_FittedShortDataAnalysis_4D() {
        if (file.IsOpen()) {
            file.Close();
        }
        if (file_hists.IsOpen()) {
            file_hists.Close();
        }
    }

    bool cutFunction() {
        return true;
    }

    void analysisEvent() {
        TH1F* missingMass = dynamic_cast<TH1F*>(file_hists.Get(nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi).c_str()));
        TCanvas* canvas = new TCanvas("canvas", "Histogram", 1200, 900);

        std::string strSubFunc = "[4]*x*x + [5]*x + [6]";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2])^2) + " + strSubFunc;

        TF1* fitFuncs = new TF1("FitFunc", strFitFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));
        TF1* subFuncs = new TF1("SubFunc", strSubFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));

        fitFuncs->SetParameter(0, ampGauss);
        fitFuncs->SetParameter(1, meanGauss);
        fitFuncs->SetParameter(2, stdDevGauss);
        fitFuncs->SetParameter(4, paramA);
        fitFuncs->SetParameter(5, paramB);
        fitFuncs->SetParameter(6, paramC);

        subFuncs->SetParameter(4, paramA);
        subFuncs->SetParameter(5, paramB);
        subFuncs->SetParameter(6, paramC);

        canvas->Divide(1, 1);
        canvas->cd(1);

        missingMass->Draw();
        fitFuncs->Draw("same");
        subFuncs->SetLineColor(kBlue);
        subFuncs->Draw("same");

        canvas->Update();
        std::string figNameCanvasMM = IMG_BUFF + "bigShort_2504/" + nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi) + ".png";
        canvas->SaveAs(figNameCanvasMM.c_str());

        delete canvas;
        delete fitFuncs;
        delete subFuncs;
    }

    void results() {}

    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            if (cutFunction()) {
                analysisEvent();
            }
        }
        results();
    }

protected:
    TFile file;
    TTree* tree;

    TFile file_hists;

    double ampGauss, meanGauss, stdDevGauss, paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, errA, errB, errC;
    int index_w, index_q2, index_cos_theta, index_phi;

    double getUpEdge(int i, int j, int k, int l) {
        const int n = 36;
        const double MM_UP_EDGES[n] = { 1.00, 1.00, 1.00, 1.02, 1.04, 1.07, 1.09, 1.10, 1.11, 1.12, 1.12, 1.10, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        return MM_UP_EDGES[(int) (j * n / NUMBER_W)];
    }

    double getDownEdge(int i, int j, int k, int l) {
        return 0.83;
    }
};
*/

/**
 * Базовый класс для работы с деревом 4д данных
*/
class MM_RootSortedDataAnalysis_4D {
public:
    MM_RootSortedDataAnalysis_4D() : file("sorted_output.root", "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << "sorted_output.root" << std::endl;
            return;
        }

        trees.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            trees[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                trees[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    trees[i][j][k].resize(NUMBER_PHI, nullptr);
                }
            }
        }
    }

    virtual ~MM_RootSortedDataAnalysis_4D() {
        if (file.IsOpen()) {
            file.Close();
        }
    }

    virtual void analysisEvent(int i, int j, int k, int l) = 0;
    virtual bool cutFunction(int i, int j, int k, int l) = 0;
    virtual void results() = 0;

    void analysisCycle() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        trees[i][j][k][l] = dynamic_cast<TTree*>(file.Get(nameHist4D_MM(i, j, k, l).c_str()));
                        if (trees[i][j][k][l]) {
                            analysisEvent(i, j, k, l);
                        }
                    }
                }
            }
        }
        results();
    }

protected:
    TFile file;
    vector<vector<vector<vector<TTree*>>>> trees;
};

/**
 * Пробегается по дереву 4д данных
 * реализует гистограммы с переменных количеством бинов
 * записывает гистограммы в файл
 * Необходимо в будущем вынести алгоритм определения оптимального количества бинов
*/
class MM_FullSortedDataAnalysis_4D : public MM_RootSortedDataAnalysis_4D {
public:
    MM_FullSortedDataAnalysis_4D() : outputFile("missingMasses_var50.root", "RECREATE") {

    };

    ~MM_FullSortedDataAnalysis_4D() {
        if (outputFile.IsOpen()) {
            outputFile.Close();
        }
    };

    void analysisEvent(int i, int j, int k, int l) override {
        trees[i][j][k][l]->SetBranchAddress("MM", &MM);
        Long64_t nEntries = trees[i][j][k][l]->GetEntries();

        std::vector<double> mm_values;
        mm_values.reserve(nEntries);  // резервируем память

        for (Long64_t f = 0; f < nEntries; f++) {
            trees[i][j][k][l]->GetEntry(f);
            if (MM >= MM_MIN && MM <= MM_MAX) {
                mm_values.push_back(MM);
            }
        }

        if (mm_values.empty()) return;

        // --------- Новый метод выбора количества бинов ---------
        size_t N = mm_values.size();

        const int B_MIN = 20;       // минимальное число бинов
        const int B_MAX = 80;      // максимальное число бинов
        const double N0 = 15000.0;  // масштаб насыщения

        double growthFactor = 1.0 - std::exp(-static_cast<double>(N) / N0);
        int optimalBins = static_cast<int>(B_MIN + (B_MAX - B_MIN) * growthFactor);
        optimalBins = std::max(B_MIN, optimalBins);  // на всякий случай

        // --------- Создание и заполнение гистограммы ---------
        TCanvas* canvas = new TCanvas("canvas", "Histogram", 1200, 900);
        std::string title = "MM: Q2 = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) +
                            " GeV^{2}, W = " + to_string_with_precision(W_MIN + STEP_W * j) + "-" + to_string_with_precision(W_MIN + STEP_W * (j + 1)) +
                            " GeV, Cos Theta: " + std::to_string(k) + ", Phi: " + std::to_string(l) + "; MM, GeV";
        std::string name = nameHist4D_MM(i, j, k, l);
        TH1F* histMM = new TH1F(name.c_str(), title.c_str(), optimalBins, MM_MIN, MM_MAX);

        for (double val : mm_values) {
            histMM->Fill(val);
        }

        outputFile.cd();
        histMM->Write();

        delete histMM;
        delete canvas;
    }

    bool cutFunction(int i, int j, int k, int l) override {
        return true;
    }

    void results() override {

    }

private:
    TFile outputFile;
    double MM;
};

/**
 * Новая реализация метода фитирования
 * В будущем разработать базовый класс
*/
class MM_TwoGaussPeakFitV2 {
public:
    MM_TwoGaussPeakFitV2(TF1* fullFitFunc) : fullFitFunc(fullFitFunc) {
        fitFunc = new TF1("func_A", strFitFunc.c_str());
    };

    ~MM_TwoGaussPeakFitV2() {
        if (fitFunc) delete fitFunc;
    }

    void setShortEdges(double downEdge, double upEdge) {
        fitFunc->SetRange(downEdge, upEdge);
    }

    void setFullEdges(double downEdge, double upEdge) {
        fullFitFunc->SetRange(downEdge, upEdge);
    }

    void setDeltaPeak(double dp) {
        deltaPeak = dp;
    }

    std::string getStrFitFunc() {
        return strFitFunc;
    }

    void fitting(TH1F* missingMass) {
        int bin_min = missingMass->FindBin(0.86);
        int bin_max = missingMass->FindBin(1.00);

        int maxBin = bin_min;
        double maxValue = missingMass->GetBinContent(bin_min);

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            double value = missingMass->GetBinContent(bin);
            if (value > maxValue) {
                maxValue = value;
                maxBin = bin;
            }
        }

        double maxPosition = missingMass->GetBinCenter(maxBin);

        fitFunc->SetParameter(0, maxValue);
        fitFunc->SetParameter(1, maxPosition);
        fitFunc->SetParameter(2, 0.03);
        fitFunc->SetParameter(3, 1.00);
        fitFunc->SetParameter(4, 0.00);
        fitFunc->SetParameter(5, 0.80);
        fitFunc->SetParameter(6, 1.60);

        fitFunc->SetParLimits(0, 0.05 * maxValue, 1.1 * maxValue);
        fitFunc->SetParLimits(1, 0.86, 1.00);
        fitFunc->SetParLimits(2, 0.0085, 0.085);
        fitFunc->SetParLimits(3, 0.33, 3.00);
        fitFunc->SetParLimits(5, 0.70, 0.95);
        fitFunc->SetParLimits(6, 0.70, 4.00);
        
        missingMass->Fit(fitFunc, "RQ");

        for (int f = 0; f < 7; f++) {
            fullFitFunc->SetParameter(f, fitFunc->GetParameter(f));
        }
        fullFitFunc->SetParameter(7, 1.00);
        fullFitFunc->SetParameter(8, deltaPeak);
        fullFitFunc->SetParameter(9, 0.03);

        fullFitFunc->SetParLimits(0, 0.7 * fitFunc->GetParameter(0), 1.3 * fitFunc->GetParameter(0));
        fullFitFunc->SetParLimits(1, 0.9 * fitFunc->GetParameter(1), 1.1 * fitFunc->GetParameter(1));
        fullFitFunc->SetParLimits(2, 0.9 * fitFunc->GetParameter(2), 1.1 * fitFunc->GetParameter(2));
        fullFitFunc->SetParLimits(3, 0.9 * fitFunc->GetParameter(3), 1.1 * fitFunc->GetParameter(3));
        fullFitFunc->SetParLimits(7, 0.00, missingMass->GetBinContent(missingMass->FindBin(deltaPeak)));
        fullFitFunc->SetParLimits(8, 0.8 * deltaPeak, 1.2 * deltaPeak);
        fullFitFunc->SetParLimits(9, 0.01, 0.1);                        

        missingMass->Fit(fullFitFunc, "RQ");
    }

    void getPrimaryParams(TF1* primaryFitFunc) {
        if (!fitFunc || !primaryFitFunc) return;

        primaryFitFunc->SetRange(fitFunc->GetXmin(), fitFunc->GetXmax());

        for (int i = 0; i < fitFunc->GetNpar(); ++i) {
            primaryFitFunc->SetParameter(i, fitFunc->GetParameter(i));
            primaryFitFunc->SetParError(i, fitFunc->GetParError(i));
        }
    }
private:
    std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + [4]*(x - [5])*(x - [6])";

    TF1* fitFunc;
    TF1* fullFitFunc;

    double deltaPeak;
};

/**
 * Работает с 4д гистограммами у которых переменное колиество бинов 
 * Проводит аппроксимацию данных второй версией фита
 * Создаёт дерево с параметрами аппроксимации и записывает его в файл
*/
class MM_FittingSortedDataAnalysis_4D : public MM_RootDataAnalysis_4D {
public:
    MM_FittingSortedDataAnalysis_4D(const char* outFileName): MM_RootDataAnalysis_4D("missingMasses_var50.root"), file(outFileName, "RECREATE"), tree("FitData", "Tree with peaks' parameters") {
        fullFitFuncs.resize(NUMBER_Q2);
        primaryFitFuncs.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            fullFitFuncs[i].resize(NUMBER_W);
            primaryFitFuncs[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                fullFitFuncs[i][j].resize(NUMBER_COS_THETA);
                primaryFitFuncs[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    fullFitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    primaryFitFuncs[i][j][k].resize(NUMBER_PHI, nullptr);
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
                        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;
                        std::string fullFuncStr = strFitFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";
                        fullFitFuncs[i][j][k][l] = new TF1(nameFullFitFunc(i, j, k, l).c_str(), fullFuncStr.c_str(), getDownEdge(i, j, k, l), getUpEdge(i, j, k, l));
                        primaryFitFuncs[i][j][k][l] = new TF1(namePrimaryFitFunc(i, j, k, l).c_str(), strFitFunc.c_str(), getDownEdge(i, j, k, l), getUpShortEdge(i, j, k, l));
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
        tree.Branch("primaryAmpGauss", &primaryAmpGauss, "primaryAmpGauss/D");
        tree.Branch("primaryMeanGauss", &primaryMeanGauss, "primaryMeanGauss/D");
        tree.Branch("primaryStdDevGauss", &primaryStdDevGauss, "primaryStdDevGauss/D");
        tree.Branch("primaryAsymGauss", &primaryAsymGauss, "primaryAsymGauss/D");
        tree.Branch("paramA", &paramA, "paramA/D");
        tree.Branch("paramB", &paramB, "paramB/D");
        tree.Branch("paramC", &paramC, "paramC/D");
        tree.Branch("primaryParamA", &primaryParamA, "primaryParamA/D");
        tree.Branch("primaryParamB", &primaryParamB, "primaryParamB/D");
        tree.Branch("primaryParamC", &primaryParamC, "primaryParamC/D");
        tree.Branch("primaryInt", &primaryInt, "primaryInt/D");
        tree.Branch("eAmpGauss", &eAmpGauss, "eAmpGauss/D");
        tree.Branch("eMeanGauss", &eMeanGauss, "eMeanGauss/D");
        tree.Branch("eStdDevGauss", &eStdDevGauss, "eStdDevGauss/D");
        tree.Branch("eAsymGauss", &eAsymGauss, "eAsymGauss/D");
        tree.Branch("eAmpGaussDelta", &eAmpGaussDelta, "eAmpGaussDelta/D");
        tree.Branch("eMeanGaussDelta", &eMeanGaussDelta, "eMeanGaussDelta/D");
        tree.Branch("eStdDevGaussDelta", &eStdDevGaussDelta, "eStdDevGaussDelta/D");
        tree.Branch("ePrimaryAmpGauss", &ePrimaryAmpGauss, "ePrimaryAmpGauss/D");
        tree.Branch("ePrimaryMeanGauss", &ePrimaryMeanGauss, "ePrimaryMeanGauss/D");
        tree.Branch("ePrimaryStdDevGauss", &ePrimaryStdDevGauss, "ePrimaryStdDevGauss/D");
        tree.Branch("ePrimaryAsymGauss", &ePrimaryAsymGauss, "ePrimaryAsymGauss/D");
        tree.Branch("errA", &errA, "errA/D");
        tree.Branch("errB", &errB, "errB/D");
        tree.Branch("errC", &errC, "errC/D");
        tree.Branch("ePrimaryParamA", &ePrimaryParamA, "ePrimaryParamA/D");
        tree.Branch("ePrimaryParamB", &ePrimaryParamB, "ePrimaryParamB/D");
        tree.Branch("ePrimaryParamC", &ePrimaryParamC, "ePrimaryParamC/D");
        tree.Branch("ePrimaryInt", &ePrimaryInt, "ePrimaryInt/D");
        tree.Branch("width_mm", &width_mm, "width_mm/D");
        tree.Branch("index_q2", &index_q2, "index_q2/I");
        tree.Branch("index_w", &index_w, "index_w/I");
        tree.Branch("index_cos_theta", &index_cos_theta, "index_cos_theta/I");
        tree.Branch("index_phi", &index_phi, "index_phi/I");
    }

    ~MM_FittingSortedDataAnalysis_4D() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        if (fullFitFuncs[i][j][k][l]) delete fullFitFuncs[i][j][k][l];
                        if (primaryFitFuncs[i][j][k][l]) delete primaryFitFuncs[i][j][k][l];
                    }
                }
            }
        }
    }

    bool cutFunction(int i, int j, int k, int l) override {
        return !(!missingMasses[i][j][k][l] || missingMasses[i][j][k][l]->GetEntries() < 1);
    }

    void analysisEvent(int i, int j, int k, int l) override {   
        MM_TwoGaussPeakFitV2 fittingProcedure(fullFitFuncs[i][j][k][l]);
        fittingProcedure.setShortEdges(getDownEdge(i, j, k, l), getUpShortEdge(i, j, k, l));
        fittingProcedure.setDeltaPeak(getDeltaPeak(i, j, k, l));
        fittingProcedure.fitting(missingMasses[i][j][k][l]);
        fittingProcedure.getPrimaryParams(primaryFitFuncs[i][j][k][l]);

        ampGauss = fullFitFuncs[i][j][k][l]->GetParameter(0);
        meanGauss = fullFitFuncs[i][j][k][l]->GetParameter(1);
        stdDevGauss = fullFitFuncs[i][j][k][l]->GetParameter(2);
        asymGauss = fullFitFuncs[i][j][k][l]->GetParameter(3);
        ampGaussDelta = fullFitFuncs[i][j][k][l]->GetParameter(7);
        meanGaussDelta = fullFitFuncs[i][j][k][l]->GetParameter(8);
        stdDevGaussDelta = fullFitFuncs[i][j][k][l]->GetParameter(9);
        primaryAmpGauss = primaryFitFuncs[i][j][k][l]->GetParameter(0);
        primaryMeanGauss = primaryFitFuncs[i][j][k][l]->GetParameter(1);
        primaryStdDevGauss = primaryFitFuncs[i][j][k][l]->GetParameter(2);
        primaryAsymGauss = primaryFitFuncs[i][j][k][l]->GetParameter(3);
        paramA = fullFitFuncs[i][j][k][l]->GetParameter(4);
        paramB = fullFitFuncs[i][j][k][l]->GetParameter(5);
        paramC = fullFitFuncs[i][j][k][l]->GetParameter(6);
        primaryParamA = primaryFitFuncs[i][j][k][l]->GetParameter(4);
        primaryParamB = primaryFitFuncs[i][j][k][l]->GetParameter(5);
        primaryParamC = primaryFitFuncs[i][j][k][l]->GetParameter(6);
        primaryInt = getPrimaryInt(missingMasses[i][j][k][l]);
        eAmpGauss = fullFitFuncs[i][j][k][l]->GetParError(0);
        eMeanGauss = fullFitFuncs[i][j][k][l]->GetParError(1);
        eStdDevGauss = fullFitFuncs[i][j][k][l]->GetParError(2);
        eAsymGauss = fullFitFuncs[i][j][k][l]->GetParError(3);
        eAmpGaussDelta = fullFitFuncs[i][j][k][l]->GetParError(7);
        eMeanGaussDelta = fullFitFuncs[i][j][k][l]->GetParError(8);
        eStdDevGaussDelta = fullFitFuncs[i][j][k][l]->GetParError(9);
        ePrimaryAmpGauss = primaryFitFuncs[i][j][k][l]->GetParError(0);
        ePrimaryMeanGauss = primaryFitFuncs[i][j][k][l]->GetParError(1);
        ePrimaryStdDevGauss = primaryFitFuncs[i][j][k][l]->GetParError(2);
        ePrimaryAsymGauss = primaryFitFuncs[i][j][k][l]->GetParError(3);
        errA = fullFitFuncs[i][j][k][l]->GetParError(4);
        errB = fullFitFuncs[i][j][k][l]->GetParError(5);
        errC = fullFitFuncs[i][j][k][l]->GetParError(6);
        ePrimaryParamA = primaryFitFuncs[i][j][k][l]->GetParError(4);
        ePrimaryParamB = primaryFitFuncs[i][j][k][l]->GetParError(5);
        ePrimaryParamC = primaryFitFuncs[i][j][k][l]->GetParError(6);
        ePrimaryInt = std::sqrt(primaryInt);
        width_mm = missingMasses[i][j][k][l]->GetBinWidth(1);
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
    vector<vector<vector<vector<TF1*>>>> fullFitFuncs;
    vector<vector<vector<vector<TF1*>>>> primaryFitFuncs;

    TFile file;
    TTree tree;

    double ampGauss, meanGauss, stdDevGauss, asymGauss, ampGaussDelta, meanGaussDelta, stdDevGaussDelta, paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss, eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta, errA, errB, errC;
    double primaryAmpGauss, primaryMeanGauss, primaryStdDevGauss, primaryInt, ePrimaryAmpGauss, ePrimaryMeanGauss, ePrimaryStdDevGauss, ePrimaryInt, width_mm, primaryAsymGauss, ePrimaryAsymGauss;
    double primaryParamA, primaryParamB, primaryParamC, ePrimaryParamA, ePrimaryParamB, ePrimaryParamC;
    int index_w, index_q2, index_cos_theta, index_phi;

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

    std::string nameFullFitFunc(int i, int j, int k, int l) {
        return "FullFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }

    std::string namePrimaryFitFunc(int i, int j, int k, int l) {
        return "PrimaryFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
    }

    double getPrimaryInt(TH1F* histMM) {
        int bin_min = histMM->FindBin(0.86);
        int bin_max = histMM->FindBin(1.00);

        double integral = 0.0;

        for (int bin = bin_min; bin <= bin_max; ++bin) {
            integral += histMM->GetBinContent(bin);
        }

        return integral;
    }
};

/**
 * По записанным в дерево параметрам аппроксимации гистограмм с переменным количеством бинов:
 * просто строит графики
*/
class MM_DrawFittedSortedDataAnalysis_4D {
public:
    MM_DrawFittedSortedDataAnalysis_4D(const char* outFileName) 
        : file(outFileName, "READ"), tree(nullptr), file_hists("missingMasses_var50.root", "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << outFileName << std::endl;
            return;
        }

        if (!file_hists.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << "missingMasses_var50.root" << std::endl;
            return;
        }

        file.GetObject("FitData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'FitData' not found in file " << outFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("ampGauss", &ampGauss);
        tree->SetBranchAddress("meanGauss", &meanGauss);
        tree->SetBranchAddress("stdDevGauss", &stdDevGauss);
        tree->SetBranchAddress("asymGauss", &asymGauss);
        tree->SetBranchAddress("ampGaussDelta", &ampGaussDelta);
        tree->SetBranchAddress("meanGaussDelta", &meanGaussDelta);
        tree->SetBranchAddress("stdDevGaussDelta", &stdDevGaussDelta);
        tree->SetBranchAddress("paramA", &paramA);
        tree->SetBranchAddress("paramB", &paramB);
        tree->SetBranchAddress("paramC", &paramC);
        tree->SetBranchAddress("index_w", &index_w);
        tree->SetBranchAddress("index_q2", &index_q2);
        tree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        tree->SetBranchAddress("index_phi", &index_phi);
    }

    ~MM_DrawFittedSortedDataAnalysis_4D() {
        if (file.IsOpen()) {
            file.Close();
        }
        if (file_hists.IsOpen()) {
            file_hists.Close();
        }
    }

    bool cutFunction() {
        return true;
    }

    void analysisEvent() {
        TH1F* missingMass = dynamic_cast<TH1F*>(file_hists.Get(nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi).c_str()));
        TCanvas* canvas = new TCanvas("canvas", "Histogram", 1200, 900);

        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc + " + [7]*exp(-0.5*((x-[8])/[9])^2)";

        TF1* fitFuncs = new TF1("FitFunc", strFitFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));
        TF1* subFuncs = new TF1("SubFunc", strSubFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpEdge(index_q2, index_w, index_cos_theta, index_phi));

        fitFuncs->SetParameter(0, ampGauss);
        fitFuncs->SetParameter(1, meanGauss);
        fitFuncs->SetParameter(2, stdDevGauss);
        fitFuncs->SetParameter(3, asymGauss);
        fitFuncs->SetParameter(4, paramA);
        fitFuncs->SetParameter(5, paramB);
        fitFuncs->SetParameter(6, paramC);
        fitFuncs->SetParameter(7, ampGaussDelta);
        fitFuncs->SetParameter(8, meanGaussDelta);
        fitFuncs->SetParameter(9, stdDevGaussDelta);

        subFuncs->SetParameter(4, paramA);
        subFuncs->SetParameter(5, paramB);
        subFuncs->SetParameter(6, paramC);

        canvas->Divide(1, 1);
        canvas->cd(1);

        missingMass->Draw();
        fitFuncs->Draw("same");
        subFuncs->SetLineColor(kBlue);
        subFuncs->Draw("same");

        canvas->Update();
        std::string figNameCanvasMM = IMG_BUFF + "bigVar50/" + nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi) + ".png";
        canvas->SaveAs(figNameCanvasMM.c_str());

        delete canvas;
        delete fitFuncs;
        delete subFuncs;
    }

    void results() {}

    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            if (cutFunction()) {
                analysisEvent();
            }
        }
        results();
    }

private:
    TFile file;
    TTree* tree;

    TFile file_hists;

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
};

/**
 * То же что предыдущий, только по первой аппроксимации, то есть один симметричный пик
*/
class MM_DrawFittedShortSortedDataAnalysis_4D {
public:
    MM_DrawFittedShortSortedDataAnalysis_4D(const char* outFileName) 
        : file(outFileName, "READ"), tree(nullptr), file_hists("missingMasses_var50.root", "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << outFileName << std::endl;
            return;
        }

        if (!file_hists.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << "missingMasses_var50.root" << std::endl;
            return;
        }

        file.GetObject("FitData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'FitData' not found in file " << outFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("primaryAmpGauss", &primaryAmpGauss);
        tree->SetBranchAddress("primaryMeanGauss", &primaryMeanGauss);
        tree->SetBranchAddress("primaryStdDevGauss", &primaryStdDevGauss);
        tree->SetBranchAddress("primaryAsymGauss", &primaryAsymGauss);
        tree->SetBranchAddress("primaryParamA", &primaryParamA);
        tree->SetBranchAddress("primaryParamB", &primaryParamB);
        tree->SetBranchAddress("primaryParamC", &primaryParamC);
        tree->SetBranchAddress("index_w", &index_w);
        tree->SetBranchAddress("index_q2", &index_q2);
        tree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        tree->SetBranchAddress("index_phi", &index_phi);
    }

    ~MM_DrawFittedShortSortedDataAnalysis_4D() {
        if (file.IsOpen()) {
            file.Close();
        }
        if (file_hists.IsOpen()) {
            file_hists.Close();
        }
    }

    bool cutFunction() {
        if (primaryAmpGauss < 0.0) return false;
        if (primaryParamA * (primaryMeanGauss - primaryParamB) * (primaryMeanGauss - primaryParamC) < 0.0) return false;
        return true;
    }

    void analysisEvent() {
        TH1F* missingMass = dynamic_cast<TH1F*>(file_hists.Get(nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi).c_str()));
        TCanvas* canvas = new TCanvas("canvas", "Histogram", 1200, 900);

        std::string strSubFunc = "[4]*(x - [5])*(x - [6])";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;

        TF1* fitFuncs = new TF1("FitFunc", strFitFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpShortEdge(index_q2, index_w, index_cos_theta, index_phi));
        TF1* subFuncs = new TF1("SubFunc", strSubFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpShortEdge(index_q2, index_w, index_cos_theta, index_phi));

        fitFuncs->SetParameter(0, primaryAmpGauss);
        fitFuncs->SetParameter(1, primaryMeanGauss);
        fitFuncs->SetParameter(2, primaryStdDevGauss);
        fitFuncs->SetParameter(3, primaryAsymGauss);
        fitFuncs->SetParameter(4, primaryParamA);
        fitFuncs->SetParameter(5, primaryParamB);
        fitFuncs->SetParameter(6, primaryParamC);

        subFuncs->SetParameter(4, primaryParamA);
        subFuncs->SetParameter(5, primaryParamB);
        subFuncs->SetParameter(6, primaryParamC);

        canvas->Divide(1, 1);
        canvas->cd(1);

        missingMass->Draw();
        fitFuncs->Draw("same");
        subFuncs->SetLineColor(kBlue);
        subFuncs->Draw("same");

        canvas->Update();
        std::string figNameCanvasMM = IMG_BUFF + "bigShortVar50/" + nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi) + ".png";
        canvas->SaveAs(figNameCanvasMM.c_str());

        delete canvas;
        delete fitFuncs;
        delete subFuncs;
    }

    void results() {}

    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            if (cutFunction()) {
                analysisEvent();
            }
        }
        results();
    }

private:
    TFile file;
    TTree* tree;

    TFile file_hists;

    double primaryAmpGauss, primaryMeanGauss, primaryStdDevGauss, primaryAsymGauss, primaryParamA, primaryParamB, primaryParamC;
    int index_w, index_q2, index_cos_theta, index_phi;

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
};

/**
 * презентация среза данных
 * не реализовано
*/
class MM_FittedSortedDataAnalysis_4D_PRESENTATION {
public:
    MM_FittedSortedDataAnalysis_4D_PRESENTATION(const char* outFileName) 
        : file(outFileName, "READ"), tree(nullptr), file_hists("missingMasses_var50.root", "READ") {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << outFileName << std::endl;
            return;
        }

        if (!file_hists.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << "missingMasses_var50.root" << std::endl;
            return;
        }

        file.GetObject("FitData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'FitData' not found in file " << outFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("ampGauss", &ampGauss);
        tree->SetBranchAddress("meanGauss", &meanGauss);
        tree->SetBranchAddress("stdDevGauss", &stdDevGauss);
        tree->SetBranchAddress("asymGauss", &asymGauss);
        tree->SetBranchAddress("ampGaussDelta", &ampGaussDelta);
        tree->SetBranchAddress("meanGaussDelta", &meanGaussDelta);
        tree->SetBranchAddress("stdDevGaussDelta", &stdDevGaussDelta);
        tree->SetBranchAddress("paramA", &paramA);
        tree->SetBranchAddress("paramB", &paramB);
        tree->SetBranchAddress("paramC", &paramC);
        /*
        tree->SetBranchAddress("primaryAmpGauss", &primaryAmpGauss);
        tree->SetBranchAddress("primaryMeanGauss", &primaryMeanGauss);
        tree->SetBranchAddress("primaryStdDevGauss", &primaryStdDevGauss);
        tree->SetBranchAddress("primaryAsymGauss", &primaryAsymGauss);
        tree->SetBranchAddress("primaryParamA", &primaryParamA);
        tree->SetBranchAddress("primaryParamB", &primaryParamB);
        tree->SetBranchAddress("primaryParamC", &primaryParamC);
        */
        tree->SetBranchAddress("eAmpGauss", &eAmpGauss);
        tree->SetBranchAddress("eMeanGauss", &eMeanGauss);
        tree->SetBranchAddress("eStdDevGauss", &eStdDevGauss);
        tree->SetBranchAddress("eAsymGauss", &eAsymGauss);
        tree->SetBranchAddress("eAmpGaussDelta", &eAmpGaussDelta);
        tree->SetBranchAddress("eMeanGaussDelta", &eMeanGaussDelta);
        tree->SetBranchAddress("eStdDevGaussDelta", &eStdDevGaussDelta);
        tree->SetBranchAddress("errA", &errA);
        tree->SetBranchAddress("errB", &errB);
        tree->SetBranchAddress("errC", &errC);
        /*
        tree->SetBranchAddress("ePrimaryAmpGauss", &ePrimaryAmpGauss);
        tree->SetBranchAddress("ePrimaryMeanGauss", &ePrimaryMeanGauss);
        tree->SetBranchAddress("ePrimaryStdDevGauss", &ePrimaryStdDevGauss);
        tree->SetBranchAddress("ePrimaryAsymGauss", &ePrimaryAsymGauss);
        tree->SetBranchAddress("ePrimaryParamA", &ePrimaryParamA);
        tree->SetBranchAddress("ePrimaryParamB", &ePrimaryParamB);
        tree->SetBranchAddress("ePrimaryParamC", &ePrimaryParamC);
        */
        tree->SetBranchAddress("index_w", &index_w);
        tree->SetBranchAddress("index_q2", &index_q2);
        tree->SetBranchAddress("index_cos_theta", &index_cos_theta);
        tree->SetBranchAddress("index_phi", &index_phi);
        /*
        tree->SetBranchAddress("primaryInt", &primaryInt);
        tree->SetBranchAddress("ePrimaryInt", &ePrimaryInt);
        tree->SetBranchAddress("width_mm", &width_mm);
        */
       std::cout << "C 1" << "\n";

        canvasHist.resize(100, nullptr);
        canvasYeild.resize(100, nullptr);
        graphs.resize(4);
        PrimaryGraphs.resize(4);
        for (auto& vec : graphs) {
            vec.resize(100);
            for (auto& g : vec) {
                g = new TGraphErrors(NUMBER_W);
            }
        }
        for (auto& vec : PrimaryGraphs) {
            vec.resize(100);
            for (auto& g : vec) {
                g = new TGraphErrors(NUMBER_W);
            }
        }
        std::cout << "C 2" << "\n";
        for (int i = 0; i < 100; i++) {
            canvasHist[i] = new TCanvas(("canvasH" + std::to_string(i + 1)).c_str(), "HistogramH", 4800, 3600);
            canvasHist[i]->Divide(4, 3);
            canvasYeild[i] = new TCanvas(("canvasY" + std::to_string(i + 1)).c_str(), "HistogramY", 1200, 900);
            canvasYeild[i]->Divide(2, 2);
            for (int j = 0; j < 4; j++) {
                std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * j], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * j + 1], 1) + " GeV^{2}, cell-" + std::to_string(i + 1) + "; W, GeV; Yield";
                std::string PrTitle = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * j], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * j + 1], 1) + " GeV^{2}, PR-cell-" + std::to_string(i + 1) + "; W, GeV; Yield";
                graphs[j][i]->SetTitle(title.c_str());
                PrimaryGraphs[j][i]->SetTitle(PrTitle.c_str());
            }
        }

        std::cout << "C DONE" << "\n";
    }

    ~MM_FittedSortedDataAnalysis_4D_PRESENTATION() {
        if (file.IsOpen()) {
            file.Close();
        }
        if (file_hists.IsOpen()) {
            file_hists.Close();
        }

        for (int i = 0; i < 100; i++) {
            if (canvasHist[i]) delete canvasHist[i];
            if (canvasYeild[i]) delete canvasYeild[i];
            for (int j = 0; j < 4; j++) {
                if (graphs[j][i]) delete graphs[j][i];
                if (PrimaryGraphs[j][i]) delete PrimaryGraphs[j][i];
            }
        }
    }

    bool cutFunction() {
        if (primaryInt < 2) return false;
        if (primaryAmpGauss < 0.0) return false;
        if (primaryParamA * (primaryMeanGauss - primaryParamB) * (primaryMeanGauss - primaryParamC) < 0.0) return false;
        return true;
    }

    void analysisEvent() {
        TH1F* missingMass = dynamic_cast<TH1F*>(file_hists.Get(nameHist4D_MM(index_q2, index_w, index_cos_theta, index_phi).c_str()));

        std::string strSubFunc = "[4]*(x - [5])*(x- [6])";
        std::string strFitFunc = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + strSubFunc;

        TF1* fitFuncs = new TF1("FitFunc", strFitFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpShortEdge(index_q2, index_w, index_cos_theta, index_phi));
        TF1* subFuncs = new TF1("SubFunc", strSubFunc.c_str(), getDownEdge(index_q2, index_w, index_cos_theta, index_phi), getUpShortEdge(index_q2, index_w, index_cos_theta, index_phi));

        fitFuncs->SetParameter(0, primaryAmpGauss);
        fitFuncs->SetParameter(1, primaryMeanGauss);
        fitFuncs->SetParameter(2, primaryStdDevGauss);
        fitFuncs->SetParameter(3, primaryAsymGauss);
        fitFuncs->SetParameter(4, primaryParamA);
        fitFuncs->SetParameter(5, primaryParamB);
        fitFuncs->SetParameter(6, primaryParamC);

        subFuncs->SetParameter(4, primaryParamA);
        subFuncs->SetParameter(5, primaryParamB);
        subFuncs->SetParameter(6, primaryParamC);

        int numberRes = (index_w == 2) ? 1 : (index_w == 7) ? 2 : (index_w == 11) ? 3 : -1;
        if (numberRes > 0) {
            canvasHist[10 * index_cos_theta + index_phi]->cd(4 * (numberRes - 1) + index_q2 + 1);
            missingMass->Draw();
            fitFuncs->Draw("same");
            subFuncs->SetLineColor(kBlue);
            subFuncs->Draw("same");
        }
        
        double yPoint = sqrt(M_PI / 2) * primaryAmpGauss * abs(primaryStdDevGauss) * (1 + abs(primaryAsymGauss)) / width_mm;
        double yErr = yPoint * sqrt(pow(ePrimaryAmpGauss / primaryAmpGauss, 2) + pow(ePrimaryStdDevGauss / primaryStdDevGauss, 2) + pow(ePrimaryAsymGauss / (1 + primaryAsymGauss), 2));

        graphs[index_q2][10 * index_cos_theta + index_phi]->SetPoint(index_w, W_MIN + (index_w + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoint);
        graphs[index_q2][10 * index_cos_theta + index_phi]->SetPointError(index_w, 0.5 * (W_MAX - W_MIN) / NUMBER_W, yErr);

        PrimaryGraphs[index_q2][10 * index_cos_theta + index_phi]->SetPoint(index_w, W_MIN + (index_w + 0.5) * (W_MAX - W_MIN) / NUMBER_W, primaryInt);
        PrimaryGraphs[index_q2][10 * index_cos_theta + index_phi]->SetPointError(index_w, 0.5 * (W_MAX - W_MIN) / NUMBER_W, ePrimaryInt);
    }

    void results() {
        for (int i = 0; i < 100; i++) {
            canvasHist[i]->Update();
            std::string figNameCanvasRes = IMG_BUFF + "presentationVar50/resonances_" + std::to_string(i + 1) + ".png";
            canvasHist[i]->SaveAs(figNameCanvasRes.c_str());

            for (int j = 0; j < 4; j++) {
                // Настройка и отрисовка основного графика (синий)
                graphs[j][i]->SetLineWidth(2);
                graphs[j][i]->SetLineStyle(1);
                graphs[j][i]->SetLineColor(kBlue);

                // Настройка и отрисовка вторичного графика (зелёный)
                PrimaryGraphs[j][i]->SetLineWidth(2);
                PrimaryGraphs[j][i]->SetLineStyle(1);
                PrimaryGraphs[j][i]->SetLineColor(kGreen + 2);  // kGreen или kGreen+2 — более насыщенный зелёный

                canvasYeild[i]->cd(j + 1);
                graphs[j][i]->Draw("APC");          // Основной график
                PrimaryGraphs[j][i]->Draw("PC SAME"); // Второй график поверх

                // (опционально: можно добавить легенду, если нужно)
            }

            canvasYeild[i]->Update();
            std::string figNameCanvasYeild = IMG_BUFF + "presentationVar50/yeild_" + std::to_string(i + 1) + ".png";
            canvasYeild[i]->SaveAs(figNameCanvasYeild.c_str());
        }
    }


    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            if (cutFunction()) {
                analysisEvent();
            }
        }
        results();
    }

private:
    TFile file;
    TTree* tree;

    TFile file_hists;

    vector<TCanvas*> canvasHist;
    vector<TCanvas*> canvasYeild;
    vector<vector<TGraphErrors*>> graphs;
    vector<vector<TGraphErrors*>> PrimaryGraphs;

    double ampGauss, meanGauss, stdDevGauss, asymGauss, ampGaussDelta, meanGaussDelta, stdDevGaussDelta, paramA, paramB, paramC;
    double eAmpGauss, eMeanGauss, eStdDevGauss, eAsymGauss, eAmpGaussDelta, eMeanGaussDelta, eStdDevGaussDelta, errA, errB, errC;
    int index_w, index_q2, index_cos_theta, index_phi;
    double primaryAmpGauss, primaryMeanGauss, primaryStdDevGauss, primaryAsymGauss, primaryParamA, primaryParamB, primaryParamC, ePrimaryAmpGauss, ePrimaryMeanGauss, ePrimaryStdDevGauss, ePrimaryAsymGauss, ePrimaryParamA, ePrimaryParamB, ePrimaryParamC, primaryInt, ePrimaryInt, width_mm;

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
};