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

#include "source/dataAnalysis.h"
#include "source/cutWrapper.h"
#include "source/constants.h"
#include "source/rg-k.h"
#include "topologyCut.h"

using ElectronCuts = CommonCuts;

std::string to_string_with_precision(double value, int precision = 3) {
    std::ostringstream out;
    out.precision(precision);
    out << std::fixed << value;
    return out.str();
}

class Neutron_PiPlus_Yield_Analysis_MM : public HipoDataAnalysis {
public:
    Neutron_PiPlus_Yield_Analysis_MM(const vector<string> &dataFileNames) : HipoDataAnalysis(dataFileNames) {

        missingMasses.resize(NUMBER_Q2, vector<TH1F*>(NUMBER_W, nullptr));
        fitFuncs.resize(NUMBER_Q2, vector<TF1*>(NUMBER_W, nullptr));
        elementQ.resize(NUMBER_Q2, nullptr);
        yieldAndW.resize(NUMBER_Q2, nullptr);

        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                std::string nameHist = "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                std::string title = "MM in cell Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + "GeV^{2}, W = " + to_string_with_precision(W_MIN + STEP_W * j) + "-" + to_string_with_precision(W_MIN + STEP_W * (j + 1)) + " GeV; MM, GeV";
                missingMasses[i][j] = new TH1F(nameHist.c_str(), title.c_str(), 80, 0.8, 1.6);
                // missingMasses[i][j] = new TH1F(nameHist.c_str(), title.c_str(), 40, 0.8, 1.6);

                std::string nameFunc = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                if (j < numW_over_1) {
                    fitFuncs[i][j] = new TF1(nameFunc.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2)", missingMasses[i][j]->GetXaxis()->GetXmin(), missingMasses[i][j]->GetXaxis()->GetXmax());
                    fitFuncs[i][j]->SetParameter(0, 50000.0);
                    // fitFuncs[i][j]->SetParameter(0, 50.0);
                    fitFuncs[i][j]->SetParameter(1, 0.95);
                    fitFuncs[i][j]->SetParameter(2, 0.032);
                } else {
                    fitFuncs[i][j] = new TF1(nameFunc.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]*x*x + [7]*x + [8]", missingMasses[i][j]->GetXaxis()->GetXmin(), missingMasses[i][j]->GetXaxis()->GetXmax());
                    fitFuncs[i][j]->SetParameter(0, 28000.0);
                    // fitFuncs[i][j]->SetParameter(0, 28.0);
                    fitFuncs[i][j]->SetParameter(1, 0.95);
                    fitFuncs[i][j]->SetParameter(2, 0.035);
                    fitFuncs[i][j]->SetParameter(3, 8600.0);
                    // fitFuncs[i][j]->SetParameter(3, 8.6);
                    fitFuncs[i][j]->SetParameter(4, 1.23);
                    fitFuncs[i][j]->SetParameter(5, 0.049);
                    fitFuncs[i][j]->SetParameter(6, 0.0);
                    fitFuncs[i][j]->SetParameter(7, 0.0);
                    fitFuncs[i][j]->SetParameter(8, 0.0);
                }
            }

            std::string nameCanvasQ = "canvas_mm" + std::to_string(i + 1);
            std::string nameCanvasY = "canvas_yield" + std::to_string(i + 1);
            elementQ[i] = new TCanvas(nameCanvasQ.c_str(), "Histograms", 7200, 2700);
            yieldAndW[i] = new TCanvas(nameCanvasY.c_str(), "Histogram", 1200, 900);
        }

        vector<Cut<const DataBanks&>*> electronCuts = {
            new IsParticleNumber(),
            new IsElectronFirst(),
            new IsCorrectStatus(),
            new IsCorrectCharge(),

            new InvariantMassCut(),
            new SFCut(),
            new ElectronMomentumCut(),
            new TimeOfFlightCut(),
            new EcalFeducialEdgesCut(),
            new ElectronZVertexCut(),
            new PionComdanimationCut(),
        };

        vector<Cut<const DataBanks&, int>*> hadronCuts = {
            new HadronMomentumCut(),
            new HadronBetaCut(),
            new HadronZVertexCut(),
        };

        vector<Cut<const DataBanks&, int>*> hadronAndElectronCuts = {
            new HadronDCFiducialCut(),
        };

        vector<CutWrapper*> standartCuts = {
            new ElectronCuts(electronCuts),
            new TopologyCut(),
            new HadronCuts(hadronCuts),
            new HadronAndElectronCuts(hadronAndElectronCuts),
        };

        this->setStandartCuts(standartCuts);
    }

    ~Neutron_PiPlus_Yield_Analysis_MM() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                delete missingMasses[i][j];
                delete fitFuncs[i][j];
            }
            if (elementQ[i]) delete elementQ[i];
            if (yieldAndW[i]) delete yieldAndW[i];          
        }
    }

    void analysisEvent(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        TLorentzVector electron = TLorentzVector();
        TLorentzVector piPlus = TLorentzVector();

        for (int i = 0; i < PART->getRows(); i++) {
            switch (PART->getInt("pid", i))
            {
            case ELECTRON:
                electron.SetXYZM(PART->getFloat("px", i), PART->getFloat("py", i), PART->getFloat("pz", i), ELECTRON_MASS);
                break;
            case PI_PLUS:
                piPlus.SetXYZM(PART->getFloat("px", i), PART->getFloat("py", i), PART->getFloat("pz", i), CHARGED_PI_MASS);
                break;
            default:
                break;
            }
        }

        TLorentzVector neutron = BEAM + PROTON_TARGET - electron - piPlus;

        double w = (BEAM + PROTON_TARGET - electron).M();
        double qq = -(BEAM - electron).M2();

        int w_bin = getWBin(w);
        int qq_bin = getQ2Bin(qq);

        if (w_bin == -1 || qq_bin == -1) {
            return;
        }

        missingMasses[qq_bin][w_bin]->Fill(neutron.M());
    }

    void results() override {
        for (int i = 0; i < NUMBER_Q2; i++) {
            double yPoints[NUMBER_W];
            double eyPoints[NUMBER_W];

            elementQ[i]->Divide(6, 3);

            for (int j = 0; j < NUMBER_W; j++) {
                elementQ[i]->cd(j + 1);

                std::string fitFuncName = "FitFunc_" + std::to_string(i + 1) + "-" + std::to_string(j + 1);
                if (j < numW_over_1) {
                    fitFuncs[i][j]->SetRange(0.80, 1.00);
                } else if (j < numW_over_2) {
                    fitFuncs[i][j]->SetRange(0.80, W_MIN + j * STEP_W - 0.15);
                } else {
                    fitFuncs[i][j]->SetRange(0.80, 1.50);
                }
                missingMasses[i][j]->Fit(fitFuncName.c_str(), "R");

                missingMasses[i][j]->Draw();
                fitFuncs[i][j]->Draw("same");
            }

            elementQ[i]->Update();
            std::string fig_MM_elements_name = "MM_elements_" + std::to_string(i + 1) + ".png";
            elementQ[i]->SaveAs(fig_MM_elements_name.c_str());

            std::string nameHist = "hist" + std::to_string(i + 1);
            std::string title = "Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}; W, GeV; Yield";
            TH1D* hist = new TH1D(nameHist.c_str(), title.c_str(), NUMBER_W, W_MIN, W_MAX);

            for (int j = 0; j < NUMBER_W; j++) {
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

private:
    const double W_MIN = 1.1;
    const double W_MAX = 2.0;

    const int NUMBER_W = 18;
    const int numW_over_1 = 6;
    const int numW_over_2 = 12;
    const int NUMBER_Q2 = 4;

    const double STEP_W = (W_MAX - W_MIN) / NUMBER_W;
    const double STEPS_Q2[8] = { 0.4, 0.6, 0.6, 1.0, 1.0, 2.0, 2.0, 3.5 };


    vector<vector<TH1F*>> missingMasses;
    vector<vector<TF1*>> fitFuncs;

    vector<TCanvas*> elementQ;
    vector<TCanvas*> yieldAndW;

    int getWBin(double w) {
        int bin = (int) ((w - W_MIN) / STEP_W);
        if (bin < 0 || bin >= NUMBER_W) {
            return -1;
        }
        return bin;
    }

    int getQ2Bin(double qq) {
        for (int i = 0; i < NUMBER_Q2; i++) {
            if (qq > STEPS_Q2[2 * i] && qq < STEPS_Q2[2 * i + 1]) {
                return i;
            }
        }
        return -1;
    }

};