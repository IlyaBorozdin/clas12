#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <cmath>

#include "source/rootDataAnalysis.h"
#include "source/constants.h"

using namespace std;

std::string to_string_with_precision(double value, int precision = 3) {
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << value;
    return out.str();
}

class Neutron_PiPlus_Root_MM : public RootDataAnalysis {
public:
    Neutron_PiPlus_Root_MM(const std::string& rootFileName) : RootDataAnalysis(rootFileName) {
        missingMasses.resize(numberQQ, vector<TH1F*>(numberW, nullptr));
        fitFuncs.resize(numberQQ, vector<TF1*>(numberW, nullptr));
        elementQ_bare.resize(numberQQ, nullptr);
        // elementQ.resize(numberQQ, nullptr);
        // yieldAndW.resize(numberQQ, nullptr);

        for (int i = 0; i < numberQQ; i++) {
            for (int j = 0; j < numberW; j++) {
                std::string nameHist = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string title = "MM in cell Q^{2} = " + to_string_with_precision(QQedge[2 * i], 1) + "-" + to_string_with_precision(QQedge[2 * i + 1], 1) + " GeV^{2}, W = " + to_string_with_precision(WMin + WStep * j) + "-" + to_string_with_precision(WMin + WStep * (j + 1)) + " GeV; MM, GeV";
                missingMasses[i][j] = new TH1F(nameHist.c_str(), title.c_str(), 80, 0.8, 1.6);
                // missingMasses[i][j]->SetTitleFont(42);
                // missingMasses[i][j]->SetTitleSize(0.15);
                // missingMasses[i][j] = new TH1F(nameHist.c_str(), title.c_str(), 40, 0.8, 1.6);
            }

            std::string nameCanvasQ_bare = "canvas_mm_bare" + std::to_string(i + 1);
            // elementQ_bare[i] = new TCanvas(nameCanvasQ_bare.c_str(), "Histograms", 7200, 2700);
            elementQ_bare[i] = new TCanvas(nameCanvasQ_bare.c_str(), "Histograms", 7200, 7200);
            // std::string nameCanvasQ = "canvas_mm" + std::to_string(i + 1);
            // elementQ[i] = new TCanvas(nameCanvasQ.c_str(), "Histograms", 7200, 2700);
            // std::string nameCanvasY = "canvas_yield" + std::to_string(i + 1);
            // yieldAndW[i] = new TCanvas(nameCanvasY.c_str(), "Histograms", 1200, 900);
        }
    }

    ~Neutron_PiPlus_Root_MM() {
        for (int i = 0; i < numberQQ; i++) {
            for (int j = 0; j < numberW; j++) {
                delete missingMasses[i][j];
                delete fitFuncs[i][j];
            }
            if (elementQ_bare[i]) delete elementQ_bare[i];    
            // if (elementQ[i]) delete elementQ[i];
            // if (yieldAndW[i]) delete yieldAndW[i];       
        }
    }

    void analysisEvent() override {
        TLorentzVector electron = TLorentzVector();
        TLorentzVector piPlus = TLorentzVector();

        electron.SetXYZM(e_px, e_py, e_pz, ELECTRON_MASS);
        piPlus.SetXYZM(pi_px, pi_py, pi_pz, CHARGED_PI_MASS);

        TLorentzVector neutron = BEAM + PROTON_TARGET - electron - piPlus;
        double w = (BEAM + PROTON_TARGET - electron).M();
        double qq = -(BEAM - electron).M2();

        int w_bin = getWBin(w);
        int qq_bin = getQQBin(qq);

        if (w_bin == -1 || qq_bin == -1) {
            return;
        }

        missingMasses[qq_bin][w_bin]->Fill(neutron.M());
    }

    void results() override {
        TFile outputFile(missingMassesFile, "RECREATE");

        for (int i = 0; i < numberQQ; i++) {

            // elementQ_bare[i]->Divide(6, 3);
            elementQ_bare[i]->Divide(6, 6);

            for (int j = 0; j < numberW; j++) {
                elementQ_bare[i]->cd(j + 1);
                missingMasses[i][j]->Draw();

                missingMasses[i][j]->Write();
            }

            elementQ_bare[i]->Update();
            // std::string fig_MM_elements_name = "MM_elements_bare_50_" + std::to_string(i + 1) + ".png";
            std::string fig_MM_elements_name = "MM_elements_bare_25_" + std::to_string(i + 1) + ".png";
            elementQ_bare[i]->SaveAs(fig_MM_elements_name.c_str());
        }

        outputFile.Close();
    }

    /*
    void fitingResults() {
        TFile inputFile(missingMassesFile, "READ");
        if (inputFile.IsZombie()) {
            std::cerr << "Error: Cannot open file " << missingMassesFile << std::endl;
            return;
        }

        std::ofstream paramsFile("fitting_parameters.txt");
        if (!paramsFile.is_open()) {
            std::cerr << "Error: Cannot open file for writing fitting parameters!" << std::endl;
            return;
        }

        for (int i = 0; i < numberQQ; i++) {
            for (int j = 0; j < numberW; j++) {
                std::string histName = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                missingMasses[i][j] = dynamic_cast<TH1F*>(inputFile.Get(histName.c_str()));
                
                if (!missingMasses[i][j]) {
                    std::cerr << "Warning: Histogram " << histName << " not found!" << std::endl;
                    continue;
                }
                
                std::string nameFunc = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                
                if (j < numW_over_1) {
                    fitFuncs[i][j] = new TF1(nameFunc.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*x + [4]", missingMasses[i][j]->GetXaxis()->GetXmin(), missingMasses[i][j]->GetXaxis()->GetXmax());
                    fitFuncs[i][j]->SetParameter(0, 50000.0);
                    // fitFuncs[i][j]->SetParameter(0, 50.0);
                    fitFuncs[i][j]->SetParameter(1, 0.95);
                    fitFuncs[i][j]->SetParameter(2, 0.032);
                    fitFuncs[i][j]->SetParameter(3, 0.0);
                    fitFuncs[i][j]->SetParameter(4, 0.0);
                } else {
                    fitFuncs[i][j] = new TF1(nameFunc.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]*x*x + [7]*x + [8]", missingMasses[i][j]->GetXaxis()->GetXmin(), missingMasses[i][j]->GetXaxis()->GetXmax());
                    fitFuncs[i][j]->SetParameter(0, 50000.0);
                    // fitFuncs[i][j]->SetParameter(0, 28.0);
                    fitFuncs[i][j]->SetParameter(1, 0.95);
                    fitFuncs[i][j]->SetParameter(2, 0.035);
                    fitFuncs[i][j]->SetParameter(3, 20000.0);
                    // fitFuncs[i][j]->SetParameter(3, 8.6);
                    fitFuncs[i][j]->SetParameter(4, 1.20);
                    fitFuncs[i][j]->SetParameter(5, 0.049);
                    fitFuncs[i][j]->SetParameter(6, 0.0);
                    fitFuncs[i][j]->SetParameter(7, 0.0);
                    fitFuncs[i][j]->SetParameter(8, 0.0);
                }

                if (j < numW_over_1) {
                    fitFuncs[i][j]->SetRange(0.80, 1.05);
                } else if (j < numW_over_2) {
                    fitFuncs[i][j]->SetRange(0.80, WMin + j * WStep - 0.15);
                } else {
                    fitFuncs[i][j]->SetRange(0.80, 1.50);
                }
                

                fitFuncs[i][j] = new TF1(nameFunc.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*x*x + [4]*x + [5]", missingMasses[i][j]->GetXaxis()->GetXmin(), missingMasses[i][j]->GetXaxis()->GetXmax());
                fitFuncs[i][j]->SetParameter(0, 50000.0);
                fitFuncs[i][j]->SetParameter(1, 0.94);
                fitFuncs[i][j]->SetParameter(2, 0.032);
                fitFuncs[i][j]->SetParameter(3, 0.0);
                fitFuncs[i][j]->SetParameter(4, 0.0);
                fitFuncs[i][j]->SetParameter(5, 0.0);

                fitFuncs[i][j]->SetRange(0.82, MM_Range[j]);

                std::string input;
                while (true) {
                    missingMasses[i][j]->Fit(nameFunc.c_str(), "R");
                    TCanvas* demonstrationCanvas = new TCanvas("demo_canvas", "Histogram", 1200, 900);

                    demonstrationCanvas->Divide(1, 1);
                    demonstrationCanvas->cd(1);

                    missingMasses[i][j]->Draw();
                    fitFuncs[i][j]->Draw("same");
                    demonstrationCanvas->Update();

                    std::cout << "Введите команду (n - для продолжения, иначе введите параметры через пробел): ";
                    if (std::cin.peek() == '\n') std::cin.ignore();
                    std::getline(std::cin, input);

                    delete demonstrationCanvas;

                    if (input == "n") {
                        // Запись параметров и погрешностей
                        paramsFile << "Fit parameters for cell (" << i + 1 << ", " << j + 1 << "):\n";
                        for (int p = 0; p < fitFuncs[i][j]->GetNpar(); ++p) {
                            double param = fitFuncs[i][j]->GetParameter(p);
                            double error = fitFuncs[i][j]->GetParError(p);
                            paramsFile << "Parameter " << p << ": " << param << " ± " << error << "\n";
                        }
                        paramsFile << "\n";
                        break;
                    }

                    // Разбираем ввод пользователя
                    std::istringstream iss(input);
                    std::vector<double> values;
                    double num;
                    while (iss >> num) {
                        values.push_back(num);
                    }

                    if (!values.empty() && values.size() <= static_cast<size_t>(fitFuncs[i][j]->GetNpar())) {
                        for (size_t k = 0; k < values.size(); ++k) {
                            fitFuncs[i][j]->SetParameter(k, values[k]);
                        }
                        std::cout << "Параметры фитирования обновлены!" << std::endl;
                    } else if (!values.empty()) {
                        std::cerr << "Ошибка: введено неверное количество параметров. Ожидается от 1 до " 
                                << fitFuncs[i][j]->GetNpar() << " значений.\n";
                    }
                }
            }
        }

        inputFile.Close();
        paramsFile.close();

        for (int i = 0; i < numberQQ; i++) {
            double yPoints[numberW];
            double eyPoints[numberW];

            elementQ[i]->Divide(6, 3);

            for (int j = 0; j < numberW; j++) {
                elementQ[i]->cd(j + 1);
                missingMasses[i][j]->Draw();
                fitFuncs[i][j]->Draw("same");
            }

            elementQ[i]->Update();
            std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + ".png";
            elementQ[i]->SaveAs(fig_MM_elements_name.c_str());

            std::string nameHist = "hist" + std::to_string(i + 1);
            std::string title = "Q^{2} = " + to_string_with_precision(QQedge[2 * i], 1) + "-" + to_string_with_precision(QQedge[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TH1D* hist = new TH1D(nameHist.c_str(), title.c_str(), numberW, WMin, WMax);

            for (int j = 0; j < numberW; j++) {
                yPoints[j] = sqrt(2 * M_PI) * fitFuncs[i][j]->GetParameter(0) * abs(fitFuncs[i][j]->GetParameter(2)) / missingMasses[i][j]->GetXaxis()->GetBinWidth(1);
                eyPoints[j] = yPoints[j] * sqrt(pow(fitFuncs[i][j]->GetParError(0) / fitFuncs[i][j]->GetParameter(0), 2) + pow(fitFuncs[i][j]->GetParError(2) / fitFuncs[i][j]->GetParameter(2), 2));

                hist->SetBinContent(j + 1, yPoints[j]);
                hist->SetBinError(j + 1, eyPoints[j]);
            }

            yieldAndW[i]->Divide(1, 1);
            yieldAndW[i]->cd(1);
            hist->Draw("E");
            yieldAndW[i]->Update();
            std::string fig_yield_name = "MM_yield" + std::to_string(i + 1) + ".png";
            yieldAndW[i]->SaveAs(fig_yield_name.c_str());

            delete hist;
        }
    }
    */

private:
    const double WMin = 1.1;
    const double WMax = 2.0;

    // const int numberW = 18;
    const int numberW = 36;
    // const int numW_over_1 = 6;
    // const int numW_over_2 = 12;
    const int numberQQ = 4;

    const double WStep = (WMax - WMin) / numberW;
    const double QQedge[8] = { 0.4, 0.6, 0.6, 1.0, 1.0, 2.0, 2.0, 4.0 };
    // const double MM_Range[18] = { 0.98, 1.01, 1.05, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 };

    // const char* missingMassesFile = "missingMasses_50.root";
    const char* missingMassesFile = "missingMasses_25E.root";

    vector<vector<TH1F*>> missingMasses;
    vector<vector<TF1*>> fitFuncs;
    vector<TCanvas*> elementQ_bare;

    // vector<TCanvas*> elementQ;
    // vector<TCanvas*> yieldAndW;

    int getWBin(double w) {
        int bin = (int) ((w - WMin) / WStep);
        if (bin < 0 || bin >= numberW) {
            return -1;
        }
        return bin;
    }

    int getQQBin(double qq) {
        for (int i = 0; i < numberQQ; i++) {
            if (qq > QQedge[2 * i] && qq < QQedge[2 * i + 1]) {
                return i;
            }
        }
        return -1;
    }
};