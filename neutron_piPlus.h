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

class Neutron_PiPlus : public HipoDataAnalysis {
public:
    Neutron_PiPlus(const vector<string> &dataFileNames) : HipoDataAnalysis(dataFileNames) {

        electron_energy = new TH1F("electron_energy", "Electron Energy - E; E, GeV", 200, 0.5, 6.8);
        pi_plus_energy = new TH1F("pi_plus_energy", "#pi^{+} Energy - E; E, GeV", 200, 0, 4);
        electron_theta = new TH1F("electron_theta", "Electron #theta; #theta", 200, 0.05, 0.5);
        pi_plus_theta = new TH1F("pi_plus_theta", "#pi^{+} #theta; #theta", 200, 0, 2.5);
        electron_phi = new TH1F("electron_phi", "Electron #phi; #phi", 200, -4, 4);
        pi_plus_phi = new TH1F("pi_plus_phi", "#pi^{+} #phi; #phi", 200, -4, 4);
        electron_vz = new TH1F("electron_vz", "Electron V_{z}; V_{z}, cm", 200, -12, 6);
        pi_plus_vz = new TH1F("pi_plus_vz", "#pi^{+} V_{z}; V_{z}, cm", 200, -10, 6);

        missing_mass = new TH1F("missing_mass", "Missing Mass MM; MM, GeV", 200, 0.6, 2.6);

        electron_energy_vs_theta = new TH2F("electron_energy_vs_theta", "Electron E vs #theta; E, Gev; #theta", 200, 1, 6.5, 200, 0.09, 0.5);
        electron_energy_vs_phi = new TH2F("electron_energy_vs_phi", "Electron E vs #phi; E, Gev; #phi", 200, 1, 6.5, 200, -3.15, 3.15);
        electron_vz_vs_energy = new TH2F("electron_vz_vs_energy", "Electron V_{z} vs E; V_{z}, cm; E, GeV", 200, -12, 6, 200, 1, 6.5);
        pi_plus_energy_vs_theta = new TH2F("pi_plus_energy_vs_theta", "#pi^{+} E vs #theta; E, Gev; #theta", 200, 0.2, 4, 200, 0.2, 2.5);
        pi_plus_energy_vs_phi = new TH2F("pi_plus_energy_vs_phi", "#pi^{+} E vs #phi; E, Gev; #phi", 200, 0.2, 4, 200, -3.15, 3.15);
        pi_plus_vz_vs_energy = new TH2F("pi_plus_vz_vs_energy", "#pi^{+} V_{z} vs E; V_{z}, cm; E, GeV", 200, -10, 6, 200, 0.2, 4);

        canvas_mm = new TCanvas("canvas_mm", "Histograms", 1200, 900);
        canvas_electron_dists = new TCanvas("canvas_electron_dists", "Histograms", 1200, 1200);
        canvas_pi_plus_dists = new TCanvas("canvas_pi_plus_dists", "Histograms", 1200, 1200);
        canvas_complex_dists = new TCanvas("canvas_complex_dists", "Histograms", 1200, 1200);

        vector<Cut<const DataBanks&>*> electronCuts = {
            new IsParticleNumber(),
            new IsElectronFirst(),
            new IsCorrectStatus(),
            new IsCorrectCharge(),

            new InvariantMassCut(),
            new SFCut(),
            new ElectronMomentumCut(),
            new TimeOfFlightCut(),
            new EcalFeducialCut(),
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

    ~Neutron_PiPlus() {
        delete electron_energy;
        delete pi_plus_energy;
        delete electron_theta;
        delete pi_plus_theta;
        delete electron_phi;
        delete pi_plus_phi;
        delete electron_vz;
        delete pi_plus_vz;

        delete missing_mass;

        delete electron_energy_vs_theta;
        delete electron_energy_vs_phi;
        delete electron_vz_vs_energy;
        delete pi_plus_energy_vs_theta;
        delete pi_plus_energy_vs_phi;
        delete pi_plus_vz_vs_energy;

        delete canvas_mm;
        delete canvas_electron_dists;
        delete canvas_pi_plus_dists;
        delete canvas_complex_dists;
    }

    void analysisEvent(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        TLorentzVector electron = TLorentzVector();
        TLorentzVector piPlus = TLorentzVector();

        double electronVz, piPlusVz;

        for (int i = 0; i < PART->getRows(); i++) {
            switch (PART->getInt("pid", i))
            {
            case ELECTRON:
                electron.SetXYZM(PART->getFloat("px", i), PART->getFloat("py", i), PART->getFloat("pz", i), ELECTRON_MASS);
                electronVz = PART->getFloat("vz", i);
                break;
            case PI_PLUS:
                piPlus.SetXYZM(PART->getFloat("px", i), PART->getFloat("py", i), PART->getFloat("pz", i), CHARGED_PI_MASS);
                piPlusVz = PART->getFloat("vz", i);
                break;
            default:
                break;
            }
        }

        TLorentzVector neutron = BEAM + PROTON_TARGET - electron - piPlus;

        electron_energy->Fill(electron.E());
        pi_plus_energy->Fill(piPlus.E());
        electron_theta->Fill(electron.Theta());
        pi_plus_theta->Fill(piPlus.Theta());
        electron_phi->Fill(electron.Phi());
        pi_plus_phi->Fill(piPlus.Phi());
        electron_vz->Fill(electronVz);
        pi_plus_vz->Fill(piPlusVz);

        missing_mass->Fill(neutron.M());

        electron_energy_vs_theta->Fill(electron.E(), electron.Theta());
        electron_energy_vs_phi->Fill(electron.E(), electron.Phi());
        electron_vz_vs_energy->Fill(electronVz, electron.E());
        pi_plus_energy_vs_theta->Fill(piPlus.E(), piPlus.Theta());
        pi_plus_energy_vs_phi->Fill(piPlus.E(), piPlus.Phi());
        pi_plus_vz_vs_energy->Fill(piPlusVz, piPlus.E());
    }

    void results() override {
        TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*x*x + [4]*x + [5] + [6]*exp(-0.5*((x-[7])/[8])^2)", missing_mass->GetXaxis()->GetXmin(), missing_mass->GetXaxis()->GetXmax());
        fitFunc->SetParameter(0, 2637.0);
        fitFunc->SetParameter(1, 0.95);
        fitFunc->SetParameter(2, 0.039);
        fitFunc->SetParameter(3, -432.36);
        fitFunc->SetParameter(4, 1518.0);
        fitFunc->SetParameter(5, -701.0);
        fitFunc->SetParameter(6, 900.0);
        fitFunc->SetParameter(7, 1.2);
        fitFunc->SetParameter(8, 0.05);

        missing_mass->Fit("fitFunc", "Q");

        canvas_mm->Divide(1, 1);
        canvas_mm->cd(1);
        missing_mass->Draw();
        fitFunc->Draw("same");
        canvas_mm->Update();
        canvas_mm->SaveAs("canvas_n_pi_plus_MM.png");

        printf("Gaussian 1 Parameters:\n");
        printf("Amplitude: %f\n", fitFunc->GetParameter(0));
        printf("Mean: %f\n", fitFunc->GetParameter(1));
        printf("Mean Error: %f\n", fitFunc->GetParError(1));
        printf("Sigma: %f\n", fitFunc->GetParameter(2));

        printf("Gaussian 2 Parameters:\n");
        printf("Amplitude: %f\n", fitFunc->GetParameter(6));
        printf("Mean: %f\n", fitFunc->GetParameter(7));
        printf("Mean Error: %f\n", fitFunc->GetParError(7));
        printf("Sigma: %f\n", fitFunc->GetParameter(8));

        printf("\nQuadratic Polynomial Parameters:\n");
        printf("a: %f\n", fitFunc->GetParameter(3));
        printf("b: %f\n", fitFunc->GetParameter(4));
        printf("c: %f\n", fitFunc->GetParameter(5));

        delete fitFunc;

        canvas_electron_dists->Divide(2, 2);
        canvas_electron_dists->cd(1);
        electron_energy->Draw();
        canvas_electron_dists->cd(2);
        electron_theta->Draw();
        canvas_electron_dists->cd(3);
        electron_phi->Draw();
        canvas_electron_dists->cd(4);
        electron_vz->Draw();
        canvas_electron_dists->Update();
        canvas_electron_dists->SaveAs("canvas_electron_dists.png");

        canvas_pi_plus_dists->Divide(2, 2);
        canvas_pi_plus_dists->cd(1);
        pi_plus_energy->Draw();
        canvas_pi_plus_dists->cd(2);
        pi_plus_theta->Draw();
        canvas_pi_plus_dists->cd(3);
        pi_plus_phi->Draw();
        canvas_pi_plus_dists->cd(4);
        pi_plus_vz->Draw();
        canvas_pi_plus_dists->Update();
        canvas_pi_plus_dists->SaveAs("canvas_pi_plus_dists.png");

        canvas_complex_dists->Divide(3, 2);
        canvas_complex_dists->cd(1);
        electron_energy_vs_theta->SetStats(0);
        electron_energy_vs_theta->Draw("colz");
        canvas_complex_dists->cd(2);
        electron_energy_vs_phi->SetStats(0);
        electron_energy_vs_phi->Draw("colz");
        canvas_complex_dists->cd(3);
        electron_vz_vs_energy->SetStats(0);
        electron_vz_vs_energy->Draw("colz");
        canvas_complex_dists->cd(4);
        pi_plus_energy_vs_theta->SetStats(0);
        pi_plus_energy_vs_theta->Draw("colz");
        canvas_complex_dists->cd(5);
        pi_plus_energy_vs_phi->SetStats(0);
        pi_plus_energy_vs_phi->Draw("colz");
        canvas_complex_dists->cd(6);
        pi_plus_vz_vs_energy->SetStats(0);
        pi_plus_vz_vs_energy->Draw("colz");
        canvas_complex_dists->Update();
        canvas_complex_dists->SaveAs("canvas_complex_dists.png");
    }

private:
    TH1F* electron_energy;
    TH1F* pi_plus_energy;
    TH1F* electron_theta;
    TH1F* pi_plus_theta;
    TH1F* electron_phi;
    TH1F* pi_plus_phi;
    TH1F* electron_vz;
    TH1F* pi_plus_vz;

    TH1F* missing_mass;

    TH2F* electron_energy_vs_theta;
    TH2F* electron_energy_vs_phi;
    TH2F* electron_vz_vs_energy;
    TH2F* pi_plus_energy_vs_theta;
    TH2F* pi_plus_energy_vs_phi;
    TH2F* pi_plus_vz_vs_energy;

    TCanvas* canvas_mm;
    TCanvas* canvas_electron_dists;
    TCanvas* canvas_pi_plus_dists;
    TCanvas* canvas_complex_dists;
};


class Neutron_PiPlus_Detailed_Analysis_MM : public HipoDataAnalysis {
public:
    Neutron_PiPlus_Detailed_Analysis_MM(const vector<string> &dataFileNames) : HipoDataAnalysis(dataFileNames) {
        missingMasses.resize(numberQQ * numberW);
        fitFuncs.resize(numberQQ * numberW);

        for (int i = 0; i < numberQQ * numberW; i++) {
            std::string nameHist = "CELL " + std::to_string(i + 1);
            std::string title = "Missing Mass in CELL " + std::to_string(i + 1) + "; MM, GeV";
            missingMasses[i] = new TH1F(nameHist.c_str(), title.c_str(), 40, 0.8, 1.6);

            std::string nameFunc = "FitFunc_" + std::to_string(i + 1);
            fitFuncs[i] = new TF1(nameFunc.c_str(), "[0]*exp(-0.5*((x-[1])/[2])^2)", missingMasses[i]->GetXaxis()->GetXmin(), missingMasses[i]->GetXaxis()->GetXmax());
            fitFuncs[i]->SetParameter(0, 80.0);
            fitFuncs[i]->SetParameter(1, 0.95);
            fitFuncs[i]->SetParameter(2, 0.032);
        }

        elementQ1 = new TCanvas("canvas_mm1", "Histograms", 4800, 1800);
        elementQ2 = new TCanvas("canvas_mm2", "Histograms", 4800, 1800);
        exitAndW1 = new TCanvas("yield1", "Histogram", 1200, 900);
        exitAndW2 = new TCanvas("yield2", "Histogram", 1200, 900);

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

    ~Neutron_PiPlus_Detailed_Analysis_MM() {
        for (int i = 0; i < numberQQ * numberW; i++) {
            delete missingMasses[i];
            delete fitFuncs[i];
        }

        delete elementQ1;
        delete elementQ2;
        delete exitAndW1;
        delete exitAndW2;
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

        int cell = getCell(w, qq);

        if (cell == -1) {
            return;
        }

        missingMasses[cell]->Fill(neutron.M());
    }

    void results() override {
        double xPoints[numberW];
        double y1Points[numberW];
        double y2Points[numberW];
        double exPoints[numberW]; // Погрешности по оси x
        double ey1Points[numberW]; // Погрешности по оси y для первого графика
        double ey2Points[numberW]; // Погрешности по оси y для второго графика

        elementQ1->Divide(6, 3);
        elementQ2->Divide(6, 3);

        for (int i = 0; i < numberQQ * numberW; i++) {
            if (i < numberW) {
                elementQ1->cd(i + 1);
            }
            else {
                elementQ2->cd(i - numberW + 1);
            }

            std::string fitFuncName = "FitFunc_" + std::to_string(i + 1);
            fitFuncs[i]->SetRange(0.80, 1.00);
            missingMasses[i]->Fit(fitFuncName.c_str(), "R");

            missingMasses[i]->Draw();
            fitFuncs[i]->Draw("same");
        }

        elementQ1->Update();
        elementQ1->SaveAs("MM_elements_1.png");
        elementQ2->Update();
        elementQ2->SaveAs("MM_elements_2.png");

        for (int i = 0; i < numberW; i++) {
            xPoints[i] = WMin + (i + 0.5) * WStep;
            exPoints[i] = 0.5 * WStep;
            y1Points[i] = sqrt(2 * M_PI) * fitFuncs[i]->GetParameter(0) * abs(fitFuncs[i]->GetParameter(2)) / missingMasses[i]->GetXaxis()->GetBinWidth(1);
            y2Points[i] = sqrt(2 * M_PI) * fitFuncs[i + numberW]->GetParameter(0) * abs(fitFuncs[i + numberW]->GetParameter(2)) / missingMasses[i + numberW]->GetXaxis()->GetBinWidth(1);

            ey1Points[i] = y1Points[i] * sqrt(pow(fitFuncs[i]->GetParError(0) / fitFuncs[i]->GetParameter(0), 2) + pow(fitFuncs[i]->GetParError(2) / fitFuncs[i]->GetParameter(2), 2));
            ey2Points[i] = y1Points[i] * sqrt(pow(fitFuncs[i + numberW]->GetParError(0) / fitFuncs[i + numberW]->GetParameter(0), 2) + pow(fitFuncs[i + numberW]->GetParError(2) / fitFuncs[i + numberW]->GetParameter(2), 2));
        }

        TGraphErrors* graph1 = new TGraphErrors(numberW, xPoints, y1Points, exPoints, ey1Points);
        TGraphErrors* graph2 = new TGraphErrors(numberW, xPoints, y2Points, exPoints, ey2Points);

        graph1->SetMarkerStyle(20);
        graph1->SetTitle("Q^{2} = 0.4 - 0.6 GeV^{2}");
        graph1->GetXaxis()->SetTitle("W, Gev");
        graph1->GetYaxis()->SetTitle("Yield");

        graph2->SetMarkerStyle(20);
        graph2->SetTitle("Q^{2} = 0.6 - 1.0 GeV^{2}");
        graph2->GetXaxis()->SetTitle("W, Gev");
        graph2->GetYaxis()->SetTitle("Yield");

        exitAndW1->Divide(1, 1);
        exitAndW1->cd(1);
        graph1->Draw("AP");
        exitAndW1->Update();
        exitAndW1->SaveAs("MM_exit1.png");

        exitAndW2->Divide(1, 1);
        exitAndW2->cd(1);
        graph2->Draw("AP");
        exitAndW2->Update();
        exitAndW2->SaveAs("MM_exit2.png");

        delete graph1;
        delete graph2;
    }

private:
    const double WStep = 0.050;
    const double WMin = 1.1;
    const double WMax = 2.0;

    const double QQ1min = 0.4;
    const double QQ1max = 0.6;
    const double QQ2min = 0.6;
    const double QQ2max = 1.0;

    const int numberW = 18;
    const int numberQQ = 2;

    vector<TH1F*> missingMasses;
    TH1F* exitFunction;

    vector<TF1*> fitFuncs;

    TCanvas* elementQ1;
    TCanvas* elementQ2;
    TCanvas* exitAndW1;
    TCanvas* exitAndW2;

    int getCell(double W, double QQ) {
        if (
            W > WMax ||
            W < WMin ||
            QQ > QQ2max ||
            QQ < QQ1min
        ) {
            return -1;
        }
        return ((int) ((W - WMin) / WStep)) + (QQ > QQ2min ? numberW : 0);
    }
};