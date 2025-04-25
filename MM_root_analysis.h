#pragma once

#include <vector>
#include <iostream>
#include <fstream>

#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>

#include "n_piPlus_from_root.h"

using namespace std;

struct FitParameters {
    int i, j;
    double leftFit, rightFit;
    std::string funcFit;
    std::vector<double> params;
    std::vector<double> paramErrors;
};

class MM_Root_Analysis {
public:
    MM_Root_Analysis(int WBins) : NUMBER_W(WBins) {
        switch (NUMBER_W)
        {
        case 18:
            MM_FILE = "missingMasses_50.root";
            figureNum = 1;
            break;
        case 36:
            MM_FILE = "missingMasses_25.root";
            figureNum = 4;
            break;        
        default:
            break;
        }
        
        graphNum = NUMBER_W / figureNum;
        missingMasses.resize(NUMBER_Q2, vector<TH1F*>(NUMBER_W, nullptr));
        fitFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        fullFititFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        subFunc.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        // elementQ.resize(NUMBER_Q2, nullptr);
        elementQ.resize(NUMBER_Q2 * figureNum, nullptr);
        yieldAndW.resize(NUMBER_Q2, nullptr);

        for (int i = 0; i < NUMBER_Q2; i++) {
            //std::string nameCanvasQ = "canvas_mm" + std::to_string(i + 1);
            // elementQ[i] = new TCanvas(nameCanvasQ.c_str(), "Histograms", 7200, 2700);
            for (int j = 0; j < figureNum; j++) {
                std::string nameCanvasQ = "canvas_mm" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                elementQ[i * figureNum + j] = new TCanvas(nameCanvasQ.c_str(), "Histograms", 7200, 2700);
            }
            std::string nameCanvasY = "canvas_yield" + std::to_string(i + 1);
            yieldAndW[i] = new TCanvas(nameCanvasY.c_str(), "Histogram", 1200, 900);
        }
    }

    ~MM_Root_Analysis() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                if (fitFuncs[i][j]) delete fitFuncs[i][j];
                if (fullFititFuncs[i][j]) delete fullFititFuncs[i][j];
                if (subFunc[i][j]) delete subFunc[i][j];
            }
            // if (elementQ[i]) delete elementQ[i];
            for (int j = 0; j < figureNum; j++) {
                if (elementQ[i * figureNum + j]) delete elementQ[i * figureNum + j];
            }
            if (yieldAndW[i]) delete yieldAndW[i];        
        }
    }

    void handFitting(int i_start, int j_start) {
        
        std::ofstream paramsFile(fittingParamFile, std::ios::app);
        if (!paramsFile.is_open()) {
            std::cerr << "Error: Cannot open file for writing fitting parameters!" << std::endl;
            return;
        }
        

        TFile inputFile(MM_FILE.c_str(), "READ");
        if (inputFile.IsZombie()) {
            throw std::runtime_error("Error: Cannot open file " + std::string(MM_FILE));
        }

        for (int i = i_start; i < NUMBER_Q2; i++) {
            for (int j = j_start; j < NUMBER_W; j++) {
                std::string nameFunc = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string histName = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                missingMasses[i][j] = dynamic_cast<TH1F*>(inputFile.Get(histName.c_str()));

                if (!missingMasses[i][j]) {
                    std::cerr << "Warning: Histogram " << histName << " not found!" << std::endl;
                    // missingMasses[i][j] = new TH1F(histName.c_str(), histName.c_str(), 100, 0, 1); // Заглушка
                }

                std::string input;
                while (true) {
                    TCanvas* demonstrationCanvas = new TCanvas("demo_canvas", "Histogram", 1200, 900);
                    demonstrationCanvas->Divide(1, 1);
                    demonstrationCanvas->cd(1);
                    missingMasses[i][j]->Draw();
                    demonstrationCanvas->Update();

                    // Ввод полинома для фона
                    std::cout << "Введите полином для пьедестала (иначе, q - выход): ";
                    if (std::cin.peek() == '\n') std::cin.ignore();
                    std::getline(std::cin, input);

                    if (input == "q") {
                        std::cout << "Выход из программы...\n";
                        paramsFile.close();
                        delete demonstrationCanvas;
                        return;
                    }

                    // Проверка: не пустой ли ввод
                    if (input.empty()) {
                        std::cerr << "Ошибка: полином не должен быть пустым.\n";
                        delete demonstrationCanvas;
                        continue;
                    }

                    std::string funcFit = "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + " + input;
                    // std::string funcFit = "[0]*exp(-0.5*((x-[1])/[2])^2) + " + input;
                    try {
                        fitFuncs[i][j] = new TF1(nameFunc.c_str(), funcFit.c_str(),
                                                missingMasses[i][j]->GetXaxis()->GetXmin(), 
                                                missingMasses[i][j]->GetXaxis()->GetXmax());
                    } catch (...) {
                        std::cerr << "Ошибка: некорректное выражение для функции фитирования.\n";
                        delete demonstrationCanvas;
                        continue;
                    }

                    // Ввод границ фитирования
                    std::cout << "Введите границы фитинга (два числа через пробел): ";
                    if (std::cin.peek() == '\n') std::cin.ignore();
                    std::getline(std::cin, input);

                    std::istringstream iss_1(input);
                    double leftFit, rightFit;
                    if (!(iss_1 >> leftFit >> rightFit) || leftFit >= rightFit) {
                        std::cerr << "Ошибка: введите два корректных числа, где левый предел меньше правого.\n";
                        delete fitFuncs[i][j];
                        delete demonstrationCanvas;
                        continue;
                    }

                    fitFuncs[i][j]->SetRange(leftFit, rightFit);

                    // Ввод параметров
                    std::cout << "Введите параметры через пробел: ";
                    if (std::cin.peek() == '\n') std::cin.ignore();
                    std::getline(std::cin, input);

                    std::istringstream iss(input);
                    std::vector<double> values;
                    double num;
                    while (iss >> num) {
                        values.push_back(num);
                    }

                    int expectedParams = fitFuncs[i][j]->GetNpar();
                    if (values.size() == 0 || values.size() > static_cast<size_t>(expectedParams)) {
                        std::cerr << "Ошибка: введите от 1 до " << expectedParams << " параметров.\n";
                        delete fitFuncs[i][j];
                        delete demonstrationCanvas;
                        continue;
                    }

                    for (size_t k = 0; k < values.size(); ++k) {
                        fitFuncs[i][j]->SetParameter(k, values[k]);
                    }

                    fitFuncs[i][j]->Draw("same");
                    demonstrationCanvas->Update();

                    // Подтверждение фитирования
                    std::cout << "Введите команду (f - для продолжения, r - ввести заново): ";
                    if (std::cin.peek() == '\n') std::cin.ignore();
                    std::getline(std::cin, input);

                    if (input == "f") {
                        missingMasses[i][j]->Fit(nameFunc.c_str(), "R");
                        demonstrationCanvas->Clear();
                        missingMasses[i][j]->Draw();
                        fitFuncs[i][j]->Draw("same");
                        demonstrationCanvas->Update();
                    } else {
                        delete fitFuncs[i][j];
                        delete demonstrationCanvas;
                        continue;
                    }

                    // Сохранение параметров
                    std::cout << "Введите команду (n - для продолжения, r - ввести заново): ";
                    if (std::cin.peek() == '\n') std::cin.ignore();
                    std::getline(std::cin, input);

                    delete demonstrationCanvas;
                    if (input == "n") {
                        paramsFile << i << " " << j << "\n";
                        paramsFile << leftFit << " " << rightFit << "\n";
                        paramsFile << funcFit << "\n";
                        for (int p = 0; p < fitFuncs[i][j]->GetNpar(); ++p) {
                            paramsFile << fitFuncs[i][j]->GetParameter(p) << " ";
                        }
                        paramsFile << "\n";
                        for (int p = 0; p < fitFuncs[i][j]->GetNpar(); ++p) {
                            paramsFile << fitFuncs[i][j]->GetParError(p) << " ";
                        }
                        paramsFile << "\n";
                        break;
                    } else {
                        delete fitFuncs[i][j];
                        continue;
                    }
                }
            }
        }

        paramsFile.close();
        //prepareResults();
    }

    void prepareResults_50() {

        std::vector<FitParameters> fitParams = readFittingParameters();

        TFile inputFile(MM_FILE.c_str(), "READ");
        if (inputFile.IsZombie()) {
            throw std::runtime_error("Error: Cannot open file " + std::string(MM_FILE));
        }

        for (int i = 0; i < NUMBER_Q2; i++) {
            double yPoints[NUMBER_W];
            double eyPoints[NUMBER_W];

            elementQ[i]->Divide(6, 3);

            for (int j = 0; j < NUMBER_W; j++) {
                FitParameters fitParam = fitParams[i * NUMBER_W + j];
                elementQ[i]->cd(j + 1);

                std::string nameFunc = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string histName = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                missingMasses[i][j] = dynamic_cast<TH1F*>(inputFile.Get(histName.c_str()));

                if (!missingMasses[i][j]) {
                    std::cerr << "Warning: Histogram " << histName << " not found!" << std::endl;
                    missingMasses[i][j] = new TH1F(histName.c_str(), histName.c_str(), 100, 0, 1); // Заглушка
                }

                fitFuncs[i][j] = new TF1(nameFunc.c_str(), fitParam.funcFit.c_str(),
                                                missingMasses[i][j]->GetXaxis()->GetXmin(), 
                                                missingMasses[i][j]->GetXaxis()->GetXmax());
                fitFuncs[i][j]->SetRange(fitParam.leftFit, fitParam.rightFit);

                TF1* subFunc = new TF1("SunFunc", fitParam.funcFit.substr(64).c_str(), fitParam.leftFit, fitParam.rightFit);
                for (size_t k = 0; k < fitParam.params.size(); ++k) {
                    fitFuncs[i][j]->SetParameter(k, fitParam.params[k]);
                    if (k > 5) {
                        subFunc->SetParameter(k, fitParam.params[k]);
                    }
                }
                missingMasses[i][j]->Fit(nameFunc.c_str(), "R");

                missingMasses[i][j]->Draw();
                fitFuncs[i][j]->Draw("same");
                subFunc->SetLineColor(kBlue);
                subFunc->Draw("same");
            }

            elementQ[i]->Update();
            std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + ".png";
            elementQ[i]->SaveAs(fig_MM_elements_name.c_str());

            std::string nameHist = "YIELD " + std::to_string(i + 1);
            std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TH1D* hist = new TH1D(nameHist.c_str(), title.c_str(), NUMBER_W, W_MIN, W_MAX);
            TGraphErrors* graph = new TGraphErrors(NUMBER_W);

            for (int j = 0; j < NUMBER_W; j++) {
                yPoints[j] = sqrt(2 * M_PI) * fitFuncs[i][j]->GetParameter(0) * abs(fitFuncs[i][j]->GetParameter(2)) / missingMasses[i][j]->GetXaxis()->GetBinWidth(1);
                eyPoints[j] = yPoints[j] * sqrt(pow(fitFuncs[i][j]->GetParError(0) / fitFuncs[i][j]->GetParameter(0), 2) + pow(fitFuncs[i][j]->GetParError(2) / fitFuncs[i][j]->GetParameter(2), 2));

                hist->SetBinContent(j + 1, yPoints[j]);
                hist->SetBinError(j + 1, eyPoints[j]);

                graph->SetPoint(j, W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoints[j]);
                graph->SetPointError(j, 0.5 * (W_MAX - W_MIN) / NUMBER_W, eyPoints[j]);
            }

            graph->SetLineWidth(2);
            graph->SetLineStyle(1);

            yieldAndW[i]->Divide(1, 1);
            yieldAndW[i]->cd(1);
            hist->SetStats(0);
            hist->Draw("E");
            graph->SetLineColor(kBlue);
            graph->Draw("C SAME");  
            yieldAndW[i]->Update();
            std::string fig_yield_name = "MM_yield_" + std::to_string(i + 1) + ".png";
            yieldAndW[i]->SaveAs(fig_yield_name.c_str());

            delete hist;
        }
    }

    void prepareResults_25() {
        // const double MM_Range[36] = { 0.96, 0.98, 1.00, 1.02, 1.04, 1.07, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.12, 1.12, 1.14, 1.14, 1.14, 1.14, 1.14, 1.13, 1.13, 1.13, 1.12, 1.12, 1.12, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        const double MM_Range[36] = { 1.00, 1.00, 1.00, 1.02, 1.04, 1.07, 1.09, 1.10, 1.11, 1.12, 1.12, 1.10, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        const double MM_downEnge = 0.83;

        TFile inputFile(MM_FILE.c_str(), "READ");
        if (inputFile.IsZombie()) {
            throw std::runtime_error("Error: Cannot open file " + std::string(MM_FILE));
        }
        
        for (int i = 0; i < NUMBER_Q2; i++) {
            double yPoints[NUMBER_W];
            double eyPoints[NUMBER_W];

            // elementQ[i]->Divide(6, 6);

            for (int j = 0; j < NUMBER_W; j++) {

                if (j % graphNum == 0) {
                    elementQ[i * figureNum + (j / graphNum)]->Divide(3, 3);
                }

                // elementQ[i]->cd(j + 1);
                elementQ[i * figureNum + (j / graphNum)]->cd((j % graphNum) + 1);

                std::string nameFunc = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string nameSubFunc = "SubFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string histName = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                missingMasses[i][j] = dynamic_cast<TH1F*>(inputFile.Get(histName.c_str()));

                if (!missingMasses[i][j]) {
                    std::cerr << "Warning: Histogram " << histName << " not found!" << std::endl;
                    missingMasses[i][j] = new TH1F(histName.c_str(), histName.c_str(), 100, 0, 1); // Заглушка
                }

                std::string funcLine = "[4]*x*x + [5]*x + [6]";
                std::string funcStr = "[0]*exp(-0.5*((x-[1])/[2])^2) + " + funcLine;
                // std::string funcStr = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + funcLine;
                
                fitFuncs[i][j] = new TF1(nameFunc.c_str(), funcStr.c_str(), MM_downEnge, MM_Range[j]);
                subFunc[i][j] = new TF1("SunFunc", funcLine.c_str(), MM_downEnge, MM_Range[j]);

                int maxBin = missingMasses[i][j]->GetMaximumBin();
                double maxValue = missingMasses[i][j]->GetBinContent(maxBin);
                double maxPosition = missingMasses[i][j]->GetBinCenter(maxBin);

                fitFuncs[i][j]->SetParameter(0, maxValue);
                fitFuncs[i][j]->SetParameter(1, maxPosition);
                fitFuncs[i][j]->SetParameter(2, 0.03);
                // fitFuncs[i][j]->SetParameter(3, 1.00);

                fitFuncs[i][j]->SetParameter(4, 0.00);
                fitFuncs[i][j]->SetParameter(5, 0.00);
                fitFuncs[i][j]->SetParameter(6, 0.00);

                missingMasses[i][j]->Fit(nameFunc.c_str(), "R");

                subFunc[i][j]->SetParameter(4, fitFuncs[i][j]->GetParameter(4));
                subFunc[i][j]->SetParameter(5, fitFuncs[i][j]->GetParameter(5));
                subFunc[i][j]->SetParameter(6, fitFuncs[i][j]->GetParameter(6));

                missingMasses[i][j]->Draw();
                fitFuncs[i][j]->Draw("same");
                subFunc[i][j]->SetLineColor(kBlue);
                subFunc[i][j]->Draw("same");

                if (j % graphNum == graphNum - 1) {
                    elementQ[i * figureNum + (j / graphNum)]->Update();
                    std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + "-" + std::to_string(j / graphNum + 1) + ".png";
                    elementQ[i * figureNum + (j / graphNum)]->SaveAs(fig_MM_elements_name.c_str());
                }
            }

            // elementQ[i]->Update();
            // std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + ".png";
            // elementQ[i]->SaveAs(fig_MM_elements_name.c_str());

            std::string nameHist = "YIELD " + std::to_string(i + 1);
            std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TH1D* hist = new TH1D(nameHist.c_str(), title.c_str(), NUMBER_W, W_MIN, W_MAX);
            TGraphErrors* graph = new TGraphErrors(NUMBER_W);

            for (int j = 0; j < NUMBER_W; j++) {
                yPoints[j] = sqrt(2 * M_PI) * fitFuncs[i][j]->GetParameter(0) * abs(fitFuncs[i][j]->GetParameter(2)) / missingMasses[i][j]->GetXaxis()->GetBinWidth(1);
                eyPoints[j] = yPoints[j] * sqrt(pow(fitFuncs[i][j]->GetParError(0) / fitFuncs[i][j]->GetParameter(0), 2) + pow(fitFuncs[i][j]->GetParError(2) / fitFuncs[i][j]->GetParameter(2), 2));

                /*
                double ampGauss = fitFuncs[i][j]->GetParameter(0);
                double stdDevGauss = abs(fitFuncs[i][j]->GetParameter(2));
                double asymGauss = fitFuncs[i][j]->GetParameter(3);
                double eAmpGauss = fitFuncs[i][j]->GetParError(0);
                double eStdDevGauss = abs(fitFuncs[i][j]->GetParError(2));
                double eAsymGauss = fitFuncs[i][j]->GetParError(3);

                yPoints[j] = sqrt(M_PI / 2) * ampGauss * stdDevGauss * (1 + asymGauss);
                eyPoints[j] = yPoints[j] * sqrt(pow(eAmpGauss / ampGauss, 2) + pow(eStdDevGauss /stdDevGauss, 2) + pow(eAsymGauss / (1 + asymGauss), 2));
                */
                hist->SetBinContent(j + 1, yPoints[j]);
                hist->SetBinError(j + 1, eyPoints[j]);

                graph->SetPoint(j, W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoints[j]);
                graph->SetPointError(j, 0.5 * (W_MAX - W_MIN) / NUMBER_W, eyPoints[j]);
            }

            graph->SetLineWidth(2);
            graph->SetLineStyle(1);

            yieldAndW[i]->Divide(1, 1);
            yieldAndW[i]->cd(1);
            hist->SetStats(0);
            hist->Draw("E");
            graph->SetLineColor(kBlue);
            graph->Draw("C SAME");  
            yieldAndW[i]->Update();
            std::string fig_yield_name = "MM_yield_" + std::to_string(i + 1) + ".png";
            yieldAndW[i]->SaveAs(fig_yield_name.c_str());

            delete hist;
        }
    }

    void prepareResults_25_FullFit() {
        // const double MM_Range[36] = { 0.96, 0.98, 1.00, 1.02, 1.04, 1.07, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.10, 1.10, 1.10, 1.10, 1.12, 1.12, 1.14, 1.14, 1.14, 1.14, 1.14, 1.13, 1.13, 1.13, 1.12, 1.12, 1.12, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };
        const double MM_Range[36] = { 0.97, 0.98, 1.00, 1.02, 1.04, 1.07, 1.10, 1.10, 1.12, 1.14, 1.18, 1.23, 1.24, 1.24, 1.28, 1.30, 1.32, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34, 1.34 };
        const double D_Peak[36] = { 0.92, 0.92, 0.92, 0.92, 0.94, 0.94, 0.96, 0.96, 0.98, 1.08, 1.08, 1.12, 1.14, 1.16, 1.16, 1.18, 1.20, 1.20, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22, 1.22 };
        const double MM_downEnge = 0.80;

        TFile inputFile(MM_FILE.c_str(), "READ");
        if (inputFile.IsZombie()) {
            throw std::runtime_error("Error: Cannot open file " + std::string(MM_FILE));
        }
        
        for (int i = 0; i < NUMBER_Q2; i++) {
            double yPoints[NUMBER_W];
            double eyPoints[NUMBER_W];

            // elementQ[i]->Divide(6, 6);

            for (int j = 0; j < NUMBER_W; j++) {

                if (j % graphNum == 0) {
                    elementQ[i * figureNum + (j / graphNum)]->Divide(3, 3);
                }

                // elementQ[i]->cd(j + 1);
                elementQ[i * figureNum + (j / graphNum)]->cd((j % graphNum) + 1);

                std::string nameFunc = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string nameFullFunc = "FullFitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string nameSubFunc = "SubFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string histName = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                missingMasses[i][j] = dynamic_cast<TH1F*>(inputFile.Get(histName.c_str()));

                if (!missingMasses[i][j]) {
                    std::cerr << "Warning: Histogram " << histName << " not found!" << std::endl;
                    missingMasses[i][j] = new TH1F(histName.c_str(), histName.c_str(), 100, 0, 1); // Заглушка
                }

                std::string funcLine = "[4]*x*x + [5]*x + [6]";
                std::string funcStr = "[0]*exp(-0.5*((x-[1])/[2]*(x < [1] ? [3] : 1))^2) + " + funcLine;
                std::string fullFuncStr = funcStr + " + [7]*exp(-0.5*((x-[8])/[9])^2)";
                
                fitFuncs[i][j] = new TF1(nameFunc.c_str(), funcStr.c_str(), MM_downEnge, MM_Range[j]);
                fullFititFuncs[i][j] = new TF1(nameFullFunc.c_str(), fullFuncStr.c_str(), MM_downEnge, MM_Range[j]);
                subFunc[i][j] = new TF1(nameSubFunc.c_str(), funcLine.c_str(), MM_downEnge, MM_Range[j]);

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

                missingMasses[i][j]->Fit(nameFunc.c_str(), "R");

                // subFunc[i][j]->SetParameter(4, fitFuncs[i][j]->GetParameter(4));
                // subFunc[i][j]->SetParameter(5, fitFuncs[i][j]->GetParameter(5));
                // subFunc[i][j]->SetParameter(6, fitFuncs[i][j]->GetParameter(6));

                
                for (int k = 0; k < 7; k++) {
                    fullFititFuncs[i][j]->SetParameter(k, fitFuncs[i][j]->GetParameter(k));
                }
                fullFititFuncs[i][j]->SetParameter(7, 0.00);
                fullFititFuncs[i][j]->SetParameter(8, D_Peak[j]);
                fullFititFuncs[i][j]->SetParameter(9, 0.04);

                fullFititFuncs[i][j]->SetParLimits(0, 0.7 * fitFuncs[i][j]->GetParameter(0), 1.3 * fitFuncs[i][j]->GetParameter(0));
                fullFititFuncs[i][j]->SetParLimits(1, 0.8 * fitFuncs[i][j]->GetParameter(1), 1.2 * fitFuncs[i][j]->GetParameter(1));
                fullFititFuncs[i][j]->SetParLimits(2, 0.8 * fitFuncs[i][j]->GetParameter(2), 1.2 * fitFuncs[i][j]->GetParameter(2));
                fullFititFuncs[i][j]->SetParLimits(3, 0.8 * fitFuncs[i][j]->GetParameter(3), 1.2 * fitFuncs[i][j]->GetParameter(3));
                // fullFititFuncs[i][j]->SetParLimits(4, -2.0 * abs(fitFuncs[i][j]->GetParameter(4)), 2.0 * abs(fitFuncs[i][j]->GetParameter(4)));
                fullFititFuncs[i][j]->SetParLimits(7, 0.00, 0.2 * fitFuncs[i][j]->GetParameter(0));
                fullFititFuncs[i][j]->SetParLimits(8, 0.95, 1.26);
                fullFititFuncs[i][j]->SetParLimits(9, 0.01, 0.1);                        

                missingMasses[i][j]->Fit(nameFullFunc.c_str(), "R");

                for (int k = 4; k < 7; k++) {
                    subFunc[i][j]->SetParameter(k, fullFititFuncs[i][j]->GetParameter(k));
                }
                

                /*
                for (int k = 4; k < 7; k++) {
                    subFunc[i][j]->SetParameter(k, fitFuncs[i][j]->GetParameter(k));
                }
                */

                missingMasses[i][j]->Draw();
                fullFititFuncs[i][j]->Draw("same");
                //fitFuncs[i][j]->Draw("same");
                subFunc[i][j]->SetLineColor(kBlue);
                subFunc[i][j]->Draw("same");

                if (j % graphNum == graphNum - 1) {
                    elementQ[i * figureNum + (j / graphNum)]->Update();
                    std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + "-" + std::to_string(j / graphNum + 1) + ".png";
                    elementQ[i * figureNum + (j / graphNum)]->SaveAs(fig_MM_elements_name.c_str());
                }
            }

            // elementQ[i]->Update();
            // std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + ".png";
            // elementQ[i]->SaveAs(fig_MM_elements_name.c_str());

            std::string nameHist = "YIELD " + std::to_string(i + 1);
            std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TH1D* hist = new TH1D(nameHist.c_str(), title.c_str(), NUMBER_W, W_MIN, W_MAX);
            TGraphErrors* graph = new TGraphErrors(NUMBER_W);

            for (int j = 0; j < NUMBER_W; j++) {
                // yPoints[j] = sqrt(2 * M_PI) * fitFuncs[i][j]->GetParameter(0) * abs(fitFuncs[i][j]->GetParameter(2)) / missingMasses[i][j]->GetXaxis()->GetBinWidth(1);
                // eyPoints[j] = yPoints[j] * sqrt(pow(fitFuncs[i][j]->GetParError(0) / fitFuncs[i][j]->GetParameter(0), 2) + pow(fitFuncs[i][j]->GetParError(2) / fitFuncs[i][j]->GetParameter(2), 2));

                
                double ampGauss = fullFititFuncs[i][j]->GetParameter(0);
                double stdDevGauss = abs(fullFititFuncs[i][j]->GetParameter(2));
                double asymGauss = abs(fullFititFuncs[i][j]->GetParameter(3));
                double eAmpGauss = fullFititFuncs[i][j]->GetParError(0);
                double eStdDevGauss = fullFititFuncs[i][j]->GetParError(2);
                double eAsymGauss = fullFititFuncs[i][j]->GetParError(3);

                yPoints[j] = sqrt(M_PI / 2) * ampGauss * stdDevGauss * (1 + asymGauss) / missingMasses[i][j]->GetXaxis()->GetBinWidth(1);
                eyPoints[j] = yPoints[j] * sqrt(pow(eAmpGauss / ampGauss, 2) + pow(eStdDevGauss /stdDevGauss, 2) + pow(eAsymGauss / (1 + asymGauss), 2));
                
                hist->SetBinContent(j + 1, yPoints[j]);
                hist->SetBinError(j + 1, eyPoints[j]);

                graph->SetPoint(j, W_MIN + (j + 0.5) * (W_MAX - W_MIN) / NUMBER_W, yPoints[j]);
                graph->SetPointError(j, 0.5 * (W_MAX - W_MIN) / NUMBER_W, eyPoints[j]);
            }

            graph->SetLineWidth(2);
            graph->SetLineStyle(1);

            yieldAndW[i]->Divide(1, 1);
            yieldAndW[i]->cd(1);
            hist->SetStats(0);
            hist->Draw("E");
            graph->SetLineColor(kBlue);
            graph->Draw("C SAME");  
            yieldAndW[i]->Update();
            std::string fig_yield_name = "MM_yield_" + std::to_string(i + 1) + ".png";
            yieldAndW[i]->SaveAs(fig_yield_name.c_str());

            delete hist;
        }
    }

private:
    // const char* MM_FILE = "missingMasses_50.root";
    // const char* MM_FILE = "missingMasses_25.root";
    std::string MM_FILE;
    const char* fittingParamFile = "fitting_parameters.txt";
    // const char* fittingParamFile = "fitting_parameters_simple.txt";

    // const int NUMBER_W = 18;
    int NUMBER_W;
    const int NUMBER_Q2 = 4;
    const double STEPS_Q2[8] = { 0.4, 0.6, 0.6, 1.0, 1.0, 2.0, 2.0, 3.5 };
    const double W_MIN = 1.1;
    const double W_MAX = 2.0;
    int figureNum;
    int graphNum;

    vector<vector<TH1F*>> missingMasses;
    vector<vector<TF1*>> fitFuncs;
    vector<vector<TF1*>> fullFititFuncs;
    vector<vector<TF1*>> subFunc;
    
    vector<TCanvas*> elementQ;
    vector<TCanvas*> yieldAndW;

    std::vector<FitParameters> readFittingParameters() {
        std::ifstream inFile(fittingParamFile);
        if (!inFile.is_open()) {
            std::cerr << "Ошибка: не удалось открыть файл " << MM_FILE << " для чтения!" << std::endl;
            return {};
        }

        std::vector<FitParameters> fitData;
        FitParameters entry;
        std::string line;

        while (std::getline(inFile, line)) {
            std::istringstream iss1(line);
            if (!(iss1 >> entry.i >> entry.j)) {
                std::cerr << "Ошибка чтения индексов!" << std::endl;
                continue;
            }

            if (!std::getline(inFile, line)) break;
            std::istringstream iss2(line);
            if (!(iss2 >> entry.leftFit >> entry.rightFit)) {
                std::cerr << "Ошибка чтения границ фитирования!" << std::endl;
                continue;
            }

            if (!std::getline(inFile, entry.funcFit)) break; // Читаем строку с функцией

            if (!std::getline(inFile, line)) break;
            std::istringstream iss3(line);
            entry.params.clear();
            double value;
            while (iss3 >> value) {
                entry.params.push_back(value);
            }

            if (!std::getline(inFile, line)) break;
            std::istringstream iss4(line);
            entry.paramErrors.clear();
            while (iss4 >> value) {
                entry.paramErrors.push_back(value);
            }

            fitData.push_back(entry);
        }

        inFile.close();
        return fitData;
    }
};
