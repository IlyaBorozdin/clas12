#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <cmath>

#include "source/rootDataAnalysis.h"
#include "source/constants.h"
#include "MM_project_const.h"

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
        missingMasses.resize(NUMBER_Q2, vector<TH1F*>(NUMBER_W, nullptr));

        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                std::string title = "MM in cell Q^{2} = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) + " GeV^{2}, W = " + to_string_with_precision(W_MIN + STEP_W * j) + "-" + to_string_with_precision(W_MIN + STEP_W * (j + 1)) + " GeV; MM, GeV";
                missingMasses[i][j] = new TH1F(nameHistMM(i, j).c_str(), title.c_str(), NUMBER_MM, MM_MIN, MM_MAX);
            }
        }
    }

    ~Neutron_PiPlus_Root_MM() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                if (missingMasses[i][j]) delete missingMasses[i][j];
            }     
        }
    }

    void analysisEvent() override {
        TLorentzVector electron = TLorentzVector();
        TLorentzVector piPlus = TLorentzVector();

        electron.SetXYZM(e_px, e_py, e_pz, ELECTRON_MASS);
        piPlus.SetXYZM(pi_px, pi_py, pi_pz, CHARGED_PI_MASS);

        TLorentzVector neutron = BEAM + PROTON_TARGET - electron - piPlus;
        double w = (BEAM + PROTON_TARGET - electron).M();
        double q2 = -(BEAM - electron).M2();

        int w_bin = getWBin(w);
        int q2_bin = getQ2Bin(q2);

        if (w_bin == -1 || q2_bin == -1) {
            return;
        }

        missingMasses[q2_bin][w_bin]->Fill(neutron.M());
    }

    void results() override {
        TFile outputFile(MM_FILE.c_str(), "RECREATE");

        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                missingMasses[i][j]->Write();
            }
        }

        outputFile.Close();
    }

private:
    vector<vector<TH1F*>> missingMasses;

    int getWBin(double w) {
        int bin = (int) ((w - W_MIN) / STEP_W);
        if (bin < 0 || bin >= NUMBER_W) {
            return -1;
        }
        return bin;
    }

    int getQ2Bin(double q2) {
        for (int i = 0; i < NUMBER_Q2; i++) {
            if (q2 > STEPS_Q2[2 * i] && q2 < STEPS_Q2[2 * i + 1]) {
                return i;
            }
        }
        return -1;
    }
};

class Neutron_PiPlus_Root_4D_MM : public RootDataAnalysis {
public:
    Neutron_PiPlus_Root_4D_MM(const std::string& rootFileName) : RootDataAnalysis(rootFileName) {
        missingMasses.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; i++) {
            missingMasses[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; j++) {
                missingMasses[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    missingMasses[i][j][k].resize(NUMBER_PHI);
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        std::string title = "MM: Q2 = " + to_string_with_precision(STEPS_Q2[2 * i], 1) + "-" + to_string_with_precision(STEPS_Q2[2 * i + 1], 1) +
                                            " GeV^{2}, W = " + to_string_with_precision(W_MIN + STEP_W * j) + "-" + to_string_with_precision(W_MIN + STEP_W * (j + 1)) +
                                            " GeV, Cos Theta: " + std::to_string(k) + ", Phi: " + std::to_string(l) + "; MM, GeV";

                        std::string name = nameHist4D_MM(i, j, k, l);

                        missingMasses[i][j][k][l] = new TH1F(name.c_str(), title.c_str(), NUMBER_MM, MM_MIN, MM_MAX);
                    }
                }
            }
        }
    }

    ~Neutron_PiPlus_Root_4D_MM() {
        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        if (missingMasses[i][j][k][l]) delete missingMasses[i][j][k][l];
                    }
                }
            }
        }
    }

    void analysisEvent() override {
        TLorentzVector electron = TLorentzVector();
        TLorentzVector piPlus = TLorentzVector();

        electron.SetXYZM(e_px, e_py, e_pz, ELECTRON_MASS);
        piPlus.SetXYZM(pi_px, pi_py, pi_pz, CHARGED_PI_MASS);

        TLorentzVector neutron = BEAM + PROTON_TARGET - electron - piPlus;
        TLorentzVector gamma = BEAM - electron;
        double w = (BEAM + PROTON_TARGET - electron).M();
        double q2 = -gamma.M2();

        int w_bin = getWBin(w);
        int q2_bin = getQ2Bin(q2);

        if (w_bin == -1 || q2_bin == -1) {
            return;
        }

        TVector3 betaLab = TVector3();
        TLorentzVector piPlusCM = piPlus;

        betaLab.SetXYZ(0.0, 0.0, -gamma.P() / (gamma.E() + PROTON_MASS));
        piPlusCM.RotateZ(-electron.Phi());
        piPlusCM.RotateY(gamma.Theta());
        piPlusCM.Boost(betaLab);
        piPlusCM.SetXYZM(-piPlusCM.Px(), -piPlusCM.Py(), -piPlusCM.Pz(), CHARGED_PI_MASS);

        double cos_theta_pi = TMath::Cos(piPlusCM.Theta());
        double phi_pi = piPlusCM.Phi();
        
        int cos_theta_bin = getCosThetaBin(cos_theta_pi);
        int phi_bin = getPhiBin(phi_pi);

        if (cos_theta_bin == -1 || phi_bin == -1) {
            return;
        }

        missingMasses[q2_bin][w_bin][cos_theta_bin][phi_bin]->Fill(neutron.M());
    }

    void results() override {
        TFile outputFile(MM_FILE.c_str(), "RECREATE");

        for (int i = 0; i < NUMBER_Q2; i++) {
            for (int j = 0; j < NUMBER_W; j++) {
                for (int k = 0; k < NUMBER_COS_THETA; k++) {
                    for (int l = 0; l < NUMBER_PHI; l++) {
                        missingMasses[i][j][k][l]->Write();
                    }
                }
            }
        }

        outputFile.Close();
    }

private:
    vector<vector<vector<vector<TH1F*>>>> missingMasses;

    int getWBin(double w) {
        int bin = (int) ((w - W_MIN) / STEP_W);
        if (bin < 0 || bin >= NUMBER_W) {
            return -1;
        }
        return bin;
    }

    int getQ2Bin(double q2) {
        for (int i = 0; i < NUMBER_Q2; i++) {
            if (q2 > STEPS_Q2[2 * i] && q2 < STEPS_Q2[2 * i + 1]) {
                return i;
            }
        }
        return -1;
    }

    int getCosThetaBin(double cos_theta) {
        int bin = (int) ((cos_theta + 1.0) / STEP_COS_THETA);
        if (bin < 0 || bin >= NUMBER_COS_THETA) {
            return -1;
        }
        return bin;
    }

    int getPhiBin(double phi) {
        int bin = (int) ((phi + (phi < 0.0 ? 2 * M_PI : 0.0)) / STEP_PHI);
        if (bin < 0 || bin >= NUMBER_PHI) {
            return -1;
        }
        return bin;
    }
};