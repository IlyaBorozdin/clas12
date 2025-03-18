#pragma once

#include <functional>
#include <cmath>
#include <map>

#include <TH1F.h>
#include <TH2F.h>
#include "TGraph.h"
#include <TLine.h>
#include <TCanvas.h>

#include "dataAnalysis.h"
#include "constants.h"

std::pair<double, double> rotatePoint(double x, double y, int sector);
TGraph* getGraph(function<double(double)> func, double xmin, double xmax, int points);

template<typename... Args>
class Cut {
public:
    virtual ~Cut() {}
    virtual bool operator()(Args... args) = 0;
    virtual void check() const = 0;
};

class IsParticleNumber : public Cut<const DataBanks&> {
public:
    bool operator()(const DataBanks& banks) override {
        return banks.PART.getRows() > 1;
    }

    void check() const override {}
};

class IsElectronFirst : public Cut<const DataBanks&> {
public:
    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;
        return PART->getInt("pid", 0) == ELECTRON && PART->getInt("pid", 1) != ELECTRON;
    }

    void check() const override {}
};

class IsCorrectStatus : public Cut<const DataBanks&> {
public:
    bool operator()(const DataBanks& banks) override {
        int status = abs(banks.PART.getInt("status", 0));
        return status >= 2000 && status < 4000;
    }

    void check() const override {}
};

class IsCorrectCharge : public Cut<const DataBanks&> {
public:
    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;
        int charge = 0;

        for (int i = 0; i < PART->getRows(); i++) {
            charge += PART->getInt("charge", i);
        }

        return charge == 0;
    }

    void check() const override {}
};

class InvariantMassCut : public Cut<const DataBanks&> {
public:
    InvariantMassCut() {
        invarianMassHist = new TH1F("invarian_mass", "Invariant Mass W; W, GeV", 200, 0.0, 4.0);
    }

    ~InvariantMassCut() {
        delete invarianMassHist;
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        TLorentzVector electron = TLorentzVector();
        electron.SetXYZM(PART->getFloat("px", 0), PART->getFloat("py", 0), PART->getFloat("pz", 0), ELECTRON_MASS);

        double invMass = (BEAM + PROTON_TARGET - electron).M();
        
        invarianMassHist->Fill(invMass);

        return invMass < MAX;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_w", "Histograms", 600, 450);

        canvas->Divide(1, 1);
        canvas->cd(1);
        invarianMassHist->Draw();

        TLine *line = new TLine(MAX, 0, MAX, invarianMassHist->GetMaximum());
        line->SetLineColor(kRed);
        line->Draw();

        canvas->Update();
        canvas->SaveAs("CUT_Invarian_Mass.png");

        delete line;
        delete canvas;
    }
private:
    const double MAX = 3.0;
    TH1* invarianMassHist;
};

class SFCut : public Cut<const DataBanks&> {
public:
    SFCut() {
        energyVsMomentumHists.resize(SECTORS);
        for (int i = 0; i < SECTORS; i++) {
            std::string title = "E_{tot}/P Vs P for the electron in S" + std::to_string(i + 1) + "; P, GeV; E_{tot}/P";
            std::string name = "sf_sector_" + std::to_string(i + 1);
            energyVsMomentumHists[i] = new TH2F(name.c_str(), title.c_str(), 200, 1.0, 7.0, 200, 0.15, 0.32);
        }
    }

    ~SFCut() {
        for (auto& hist : energyVsMomentumHists) {
            delete hist;
        }
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* ECAL = &banks.ECAL;
        const hipo::bank* PART = &banks.PART;

        double E_pcal = .0;
        double E_in = .0;
        double E_out = .0;

        const int PRESTORM_LAYER = 1;
        const int IN_LAYER = 4;
        const int OUT_LAYER = 7;

        for (int i = 0; i < ECAL->getRows(); i++) {
            if (ECAL->getInt("pindex", i) == 0) {
                switch (ECAL->getInt("layer", i))
                {
                case PRESTORM_LAYER:
                    E_pcal = ECAL->getFloat("energy", i);
                    break;
                case IN_LAYER:
                    E_in = ECAL->getFloat("energy", i);
                    break;
                case OUT_LAYER:
                    E_out = ECAL->getFloat("energy", i);
                    break;
                default:
                    break;
                }
            }
        }

        double E_total = E_pcal + E_in + E_out;

        TLorentzVector electron = TLorentzVector();
        electron.SetXYZM(PART->getFloat("px", 0), PART->getFloat("py", 0), PART->getFloat("pz", 0), ELECTRON_MASS);
        double electronP = electron.P();

        double SF_lower, SF_upper;

        for (int i = 0; i < ECAL->getRows(); i++) {
            if (ECAL->getInt("pindex", i) == 0) {
                SF_lower = lowerSFCut(ECAL->getInt("sector", i))(electronP);
                SF_upper = upperSFCut(ECAL->getInt("sector", i))(electronP);

                energyVsMomentumHists[ECAL->getInt("sector", i) - 1]->Fill(electronP, E_total / electronP);
            }
        }

        return E_total / electronP > SF_lower && E_total /electronP < SF_upper;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_SF", "Histograms", 1200, 900);
        canvas->Divide(3, 2);

        vector<TGraph*> upperFuncs;
        vector<TGraph*> lowerFuncs;

        upperFuncs.resize(SECTORS);
        lowerFuncs.resize(SECTORS);
        
        for (int i = 0; i < SECTORS; i++) {
            canvas->cd(i + 1);
            energyVsMomentumHists[i]->Draw("colz");

            upperFuncs[i] = getGraph(upperSFCut(i + 1), 1.0, 6.0, 100);
            lowerFuncs[i] = getGraph(lowerSFCut(i + 1), 1.0, 6.0, 100);

            upperFuncs[i]->Draw("same");
            lowerFuncs[i]->Draw("same");

        }
        
        canvas->Update();
        canvas->SaveAs("CUT_SF.png");

        for (int i = 0; i < SECTORS; i++) {
            delete upperFuncs[i];
            delete lowerFuncs[i];
        }

        delete canvas;
    }
private:
    const int SECTORS = 6;
    vector<TH2F*> energyVsMomentumHists;

    function<double(double)> parabolaFunction(double p0, double p1, double p2) const {
        return [=](double x) -> double { return p0 + p1 * x / 1000 + p2 * x * x / 1000; };
    }

    function<double(double)> lowerSFCut(int sector) const {
        switch (sector) {
        case 1:
            return parabolaFunction(0.145, 21.6, -1.81);
            break;
        case 2:
            return parabolaFunction(0.134, 23.0, -1.68);
            break;
        case 3:
            return parabolaFunction(0.145, 21.1, -1.66);
            break;
        case 4:
            return parabolaFunction(0.152, 16.1, -1.27);
            break;
        case 5:
            return parabolaFunction(0.141, 21.0, -1.74);
            break;
        case 6:
            return parabolaFunction(0.141, 21.52, -1.70);
            break;    
        default:
            return [=](double p) -> double { return .0; };
            break;
        }
    }

    function<double(double)> upperSFCut(int sector) const {
        switch (sector) {
        case 1:
            return parabolaFunction(0.321, -4.58, -0.086);
            break;
        case 2:
            return parabolaFunction(0.300, 4.18, -0.650);
            break;
        case 3:
            return parabolaFunction(0.309, -0.712, -0.322);
            break;
        case 4:
            return parabolaFunction(0.322, -8.78, 0.434);
            break;
        case 5:
            return parabolaFunction(0.308, -2.01, -0.211);
            break;
        case 6:
            return parabolaFunction(0.306, -0.253, -0.328);
            break;    
        default:
            return [=](double p) -> double { return .0; };
            break;
        }
    }
};

class ElectronMomentumCut : public Cut<const DataBanks&> {
public:
    ElectronMomentumCut() {
        electronMomentumHist = new TH1F("electron_momentum", "Electron Momentum P; P, GeV", 200, 2.0, 7.0);
    }

    ~ElectronMomentumCut() {
        delete electronMomentumHist;
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        TLorentzVector electron = TLorentzVector();
        electron.SetXYZM(PART->getFloat("px", 0), PART->getFloat("py", 0), PART->getFloat("pz", 0), ELECTRON_MASS);
        double electronP = electron.P();

        electronMomentumHist->Fill(electronP);

        return electronP > LOW_P && electronP < UP_P;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_p", "Histograms", 600, 450);

        canvas->Divide(1, 1);
        canvas->cd(1);
        electronMomentumHist->Draw();

        TLine *lowLine = new TLine(LOW_P, 0, LOW_P, electronMomentumHist->GetMaximum());
        lowLine->SetLineColor(kRed);
        lowLine->Draw();

        TLine *upLine = new TLine(UP_P, 0, UP_P, electronMomentumHist->GetMaximum());
        upLine->SetLineColor(kRed);
        upLine->Draw();

        canvas->Update();
        canvas->SaveAs("CUT_Electron_Momentum.png");

        delete lowLine;
        delete upLine;
        delete canvas;
    }
private:
    const double LOW_P = 3.0;
    const double UP_P = BEAM.P();
    TH1* electronMomentumHist;
};

class TimeOfFlightCut : public Cut<const DataBanks&> {
public:
    TimeOfFlightCut() {
        timeOfFlightHist = new TH1F("time_of_flight", "Time of Flight; t, ns", 200, 17.0, 30.0);
    }

    ~TimeOfFlightCut() {
        delete timeOfFlightHist;
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;
        const hipo::bank* SCINT = &banks.SCINT;

        double startTime = PART->getFloat("vt", 0);
        vector<double> scintTimes;

        for (int i = 0; i < SCINT->getRows(); i++) {
            if (SCINT->getInt("pindex", i) == 0) {
                scintTimes.push_back(SCINT->getFloat("time", i));
            }
        }

        if (scintTimes.size() == 0) {
            return false;
        }

        double hitTime = *max_element(begin(scintTimes), end(scintTimes));
        double timeOfFlight = hitTime - startTime;

        timeOfFlightHist->Fill(timeOfFlight);

        return timeOfFlight > MIN && timeOfFlight < MAX;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_t", "Histograms", 600, 450);

        canvas->Divide(1, 1);
        canvas->cd(1);
        timeOfFlightHist->Draw();

        TLine *lowLine = new TLine(MIN, 0, MIN, timeOfFlightHist->GetMaximum());
        lowLine->SetLineColor(kRed);
        lowLine->Draw();

        TLine *upLine = new TLine(MAX, 0, MAX, timeOfFlightHist->GetMaximum());
        upLine->SetLineColor(kRed);
        upLine->Draw();

        canvas->Update();
        canvas->SaveAs("CUT_Time_of_Flight.png");

        delete lowLine;
        delete upLine;
        delete canvas;
    }
private:
    const double MIN = 21.0;
    const double MAX = 26.0;
    TH1* timeOfFlightHist;
};

class ElectronZVertexCut : public Cut<const DataBanks&> {
public:
    ElectronZVertexCut() {
        zVertexHist = new TH1F("electron_z_vertex", "Electron Z Vertex; v_{z}, cm", 200, -12.0, 4.0);
    }

    ~ElectronZVertexCut() {
        delete zVertexHist;
    }

    bool operator()(const DataBanks& banks) override {
        double vz = banks.PART.getFloat("vz", 0);

        zVertexHist->Fill(vz);

        return vz > MIN && vz < MAX;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_vz", "Histograms", 600, 450);

        canvas->Divide(1, 1);
        canvas->cd(1);
        zVertexHist->Draw();

        TLine *lowLine = new TLine(MIN, 0, MIN, zVertexHist->GetMaximum());
        lowLine->SetLineColor(kRed);
        lowLine->Draw();

        TLine *upLine = new TLine(MAX, 0, MAX, zVertexHist->GetMaximum());
        upLine->SetLineColor(kRed);
        upLine->Draw();

        canvas->Update();
        canvas->SaveAs("CUT_Electron_Z_Vertex.png");

        delete lowLine;
        delete upLine;
        delete canvas;
    }

private:
    const double MIN = -10.0;
    const double MAX = 2.0;
    TH1* zVertexHist;
};

class EcalFeducialCut : public Cut<const DataBanks&> {
public:
    EcalFeducialCut() {
        UvsVHist = new TH2F("U_vs_V", "Ecal V vs U; U, cm; V, cm", 200, 0.0, 450.0, 200, 0.0, 280.0);
        UvsWHist = new TH2F("U_vs_W", "Ecal W vs U; U, cm; W, cm", 200, 0.0, 450.0, 200, 0.0, 280.0);
    }

    ~EcalFeducialCut() {
        delete UvsVHist;
        delete UvsWHist;
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* ECAL = &banks.ECAL;

        double U_pcal = .0;
        double V_pcal = .0;
        double W_pcal = .0;

        const int PRESTORM_LAYER = 1;

        for (int i = 0; i < ECAL->getRows(); i++) {
            if (ECAL->getInt("pindex", i) == 0) {
                switch (ECAL->getInt("layer", i))
                {
                case PRESTORM_LAYER:
                    U_pcal = ECAL->getFloat("lu", i);
                    V_pcal = ECAL->getFloat("lv", i);
                    W_pcal = ECAL->getFloat("lw", i);
                    break;
                default:
                    break;
                }
            }
        }

        UvsVHist->Fill(U_pcal, V_pcal);
        UvsWHist->Fill(U_pcal, W_pcal);

        return U_pcal > lowU && U_pcal < upU && V_pcal > lowEdge && W_pcal > lowEdge;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_ecal_fecucial", "Histograms", 1200, 450);

        canvas->Divide(2, 1);
        canvas->cd(1);
        UvsVHist->Draw("colz");

        TLine *lowLine = new TLine(lowU, lowEdge, lowU, UvsVHist->GetMaximum());
        lowLine->SetLineColor(kRed);
        lowLine->Draw();

        TLine *upLine = new TLine(upU, lowEdge, upU, UvsVHist->GetMaximum());
        upLine->SetLineColor(kRed);
        upLine->Draw();

        TLine *edgeLine = new TLine(lowU, lowEdge, upU, lowEdge);
        edgeLine->SetLineColor(kRed);
        edgeLine->Draw();

        canvas->cd(2);
        UvsWHist->Draw("colz");
        lowLine->Draw();
        upLine->Draw();
        edgeLine->Draw();

        canvas->Update();
        canvas->SaveAs("CUT_ECAL_Feducial.png");

        delete lowLine;
        delete upLine;
        delete canvas;
    }
private:
    const double lowU = 40.0;
    const double upU = 400.0;
    const double lowEdge = 15.0;

    TH2F* UvsVHist;
    TH2F* UvsWHist;
};

class EcalFeducialEdgesCut : public Cut<const DataBanks&> {
public:
    EcalFeducialEdgesCut() {
        ecalFeducialSectorHists.resize(SECTORS);
        for (int i = 0; i < SECTORS; i++) {
            std::string title = "ECAL y Vs x in S" + std::to_string(i + 1) + "; x, cm; y, cm";
            std::string name = "ecal_sector_" + std::to_string(i + 1);
            ecalFeducialSectorHists[i] = new TH2F(name.c_str(), title.c_str(), 200, 40.0, 400.0, 200, -200.0, 200.);
        }
    }

    ~EcalFeducialEdgesCut() {
        for (auto& hist : ecalFeducialSectorHists) {
            delete hist;
        }
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* ECAL = &banks.ECAL;

        const int PRESTORM_LAYER = 1;

        for (int i = 0; i < ECAL->getRows(); i++) {
            if (ECAL->getInt("pindex", i) == 0 && ECAL->getInt("layer", i) == PRESTORM_LAYER) {
                double cx0 = ECAL->getFloat("x", i);
                double cy0 = ECAL->getFloat("y", i);
                int sector = ECAL->getInt("sector", i);

                std::pair<double, double> coordinate = rotatePoint(cx0, cy0, sector);

                double cx = coordinate.first;
                double cy = coordinate.second;

                double* params = ecalFeducialCutParams(sector);

                ecalFeducialSectorHists[sector - 1]->Fill(cx, cy);

                if (cx > params[P_S] && (cy > params[S_L] * (cx - params[T_L]) || cy < params[S_R] * (cx - params[T_R]))) {
                    delete params;
                    return false;
                }
                else if (cx < params[P_S] && (cy > params[Q_L] * (cx - params[R_L]) || cy < params[Q_R] * (cx - params[R_R]))) {
                    delete params;
                    return false;
                }

                delete params;
            }
        } 

        return true;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_ECAL_Feducial", "Histograms", 1200, 900);
        canvas->Divide(3, 2);

        vector<TLine*> upLines, lowLines, edgeUpLines, edgeLowLines;
        upLines.resize(SECTORS);
        lowLines.resize(SECTORS);
        edgeUpLines.resize(SECTORS);
        edgeLowLines.resize(SECTORS);

        for (int i = 0; i < SECTORS; i++) {
            canvas->cd(i + 1);
            ecalFeducialSectorHists[i]->Draw("colz");

            double* params = ecalFeducialCutParams(i + 1);

            upLines[i] = new TLine(params[P_S], params[S_L] * (params[P_S] - params[T_L]), 350, params[S_L] * (350 - params[T_L]));
            lowLines[i] = new TLine(params[P_S], params[S_R] * (params[P_S] - params[T_R]), 350, params[S_R] * (350 - params[T_R]));
            edgeUpLines[i] = new TLine((params[Q_L] * params[R_L] - params[Q_R] * params[R_R]) / (params[Q_L] - params[Q_R]), params[Q_L] * params[Q_R] * (params[R_L] - params[R_R]) / (params[Q_L] - params[Q_R]), params[P_S], params[Q_L] * (params[P_S] - params[R_L]));
            edgeLowLines[i] = new TLine((params[Q_L] * params[R_L] - params[Q_R] * params[R_R]) / (params[Q_L] - params[Q_R]), params[Q_L] * params[Q_R] * (params[R_L] - params[R_R]) / (params[Q_L] - params[Q_R]), params[P_S], params[Q_R] * (params[P_S] - params[R_R]));

            upLines[i]->SetLineColor(kRed);
            lowLines[i]->SetLineColor(kRed);
            edgeUpLines[i]->SetLineColor(kRed);
            edgeLowLines[i]->SetLineColor(kRed);

            upLines[i]->Draw();
            lowLines[i]->Draw();
            edgeUpLines[i]->Draw();
            edgeLowLines[i]->Draw();

            delete params;
        }

        canvas->Update();
        canvas->SaveAs("CUT_ECAL_Feducial_Edges.png");

        for (int i = 0; i < SECTORS; i++) {
            delete upLines[i];
            delete lowLines[i];
            delete edgeUpLines[i];
            delete edgeLowLines[i];
        }

        delete canvas;
    }
private:
    const int P_S = 0, T_L = 1, T_R = 2, S_L = 3, S_R = 4, R_L = 5, R_R = 6, Q_L = 7, Q_R = 8;
    const int SECTORS = 6;
    vector<TH2F*> ecalFeducialSectorHists;

    double* ecalFeducialCutParams(int sector) const {
        switch (sector)
        {
        case 1:
            return new double[9] {87, 58.7356, 58.7477, 0.582053, -0.591876, 64.9348, 65.424, 0.745578, -0.775022};
        case 2:
            return new double[9] {82, 62.8204, 51.2589, 0.544976, -0.562926, 64.7541, 54.6992, 0.606081, -0.633863};
        case 3:
            return new double[9] {85, 62.2296, 59.2357, 0.549788, -0.562246, 67.832, 63.6628, 0.729202, -0.678901};
        case 4:
            return new double[9] {77, 53.7756, 56.2415, 0.56899, -0.563726, 55.9324, 57.8931, 0.627239, -0.612458};
        case 5:
            return new double[9] {78, 58.2888, 60.8219, 0.56414, -0.568902, 55.9225, 56.5367, 0.503674, -0.454319};
        case 6:
            return new double[9] {82, 54.5822, 49.8914, 0.57343, -0.550729, 60.0997, 56.4641, 0.717899, -0.692481};
        default:
            return new double[0];
        }
    }
};

class PionComdanimationCut : public Cut<const DataBanks&> {
public:
    PionComdanimationCut() {
        pionComdanimationHist = new TH2F("ecin_vs_pcal", "ECIN vs PCAL SF e; E_{PCAL}/p_{e}; E_{ECin}/p_{e}", 200, 0.0, 0.25, 200, 0.0, 0.2);
    }

    ~PionComdanimationCut() {
        delete pionComdanimationHist;
    }

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* ECAL = &banks.ECAL;
        const hipo::bank* PART = &banks.PART;

        double E_pcal = .0;
        double E_in = .0;

        const int PRESTORM_LAYER = 1;
        const int IN_LAYER = 4;

        for (int i = 0; i < ECAL->getRows(); i++) {
            if (ECAL->getInt("pindex", i) == 0) {
                switch (ECAL->getInt("layer", i))
                {
                case PRESTORM_LAYER:
                    E_pcal = ECAL->getFloat("energy", i);
                    break;
                case IN_LAYER:
                    E_in = ECAL->getFloat("energy", i);
                    break;
                default:
                    break;
                }
            }
        }

        TLorentzVector electron = TLorentzVector();
        electron.SetXYZM(PART->getFloat("px", 0), PART->getFloat("py", 0), PART->getFloat("pz", 0), ELECTRON_MASS);
        double electronP = electron.P();

        pionComdanimationHist->Fill(E_pcal / electronP, E_in / electronP);

        return E_in / electronP > -0.84 * E_pcal / electronP + 0.17;
    }

    void check() const override {
        TCanvas* canvas = new TCanvas("canvas_pion_comdanimation", "Histogrtams", 600, 450);
        canvas->Divide(1, 1);

        canvas->cd(1);
        pionComdanimationHist->Draw("colz");

        TLine* line = new TLine(0.0, 0.17, 0.2024, 0.0);
        line->SetLineColor(kRed);
        line->Draw();

        canvas->Update();
        canvas->SaveAs("CUT_Pion_Comdanimation.png");
        
        delete line;
        delete canvas;
    }
private:
    TH2F* pionComdanimationHist;
};

class HadronMomentumCut : public Cut<const DataBanks&, int> {
public:
    HadronMomentumCut() {
        for (Pids hadron : CHARGED_HADRONS) {
            std::string title = "Hadron CD P_{" + std::to_string(hadron) + "}; P, GeV";
            std::string name = "hadron_CD_P_" + std::to_string(hadron);
            hadronMomentumCDHists[hadron] = new TH1F(name.c_str(), title.c_str(), 200, 0.0, 7.0);

            title = "Hadron FD P_{" + std::to_string(hadron) + "}; P, GeV";
            name = "hadron_FD_P_" + std::to_string(hadron);
            hadronMomentumFDHists[hadron] = new TH1F(name.c_str(), title.c_str(), 200, 0.0, 7.0);
        }
    }

    ~HadronMomentumCut() {
        for (auto it = hadronMomentumCDHists.begin(); it != hadronMomentumCDHists.end(); it++) {
            delete it->second;
        }

        for (auto it = hadronMomentumFDHists.begin(); it != hadronMomentumFDHists.end(); it++) {
            delete it->second;
        }
    }

    bool operator()(const DataBanks& banks, int i) override {
        const hipo::bank* PART = &banks.PART;

        int status = abs(PART->getInt("status", i));
        Pids pid = Pids(PART->getInt("pid", i));

        TLorentzVector hadron;
        hadron.SetXYZM(PART->getFloat("px", i), PART->getFloat("py", i), PART->getFloat("pz", i), partMass(pid));
        double momentum = hadron.P();

        if (status >= 2000 && status < 4000) {
            hadronMomentumFDHists[pid]->Fill(momentum);
            return momentum > MIN_FD && momentum < MAX;
        }
        hadronMomentumCDHists[pid]->Fill(momentum);
        return momentum > MIN_CD && momentum < MAX;
    }

    void check() const override {
        for (Pids hadron : CHARGED_HADRONS) {
            if (hadron == PI_PLUS) {
                TCanvas* canvas = new TCanvas("canvas_pi_plus_momentum", "Histogrtams", 1200, 450);
                canvas->Divide(2, 1);

                canvas->cd(1);
                hadronMomentumCDHists.at(hadron)->Draw();

                TLine *minLineCD = new TLine(MIN_CD, 0, MIN_CD, hadronMomentumCDHists.at(hadron)->GetMaximum());
                minLineCD->SetLineColor(kRed);
                minLineCD->Draw();

                TLine *maxLineCD = new TLine(MAX, 0, MAX, hadronMomentumCDHists.at(hadron)->GetMaximum());
                maxLineCD->SetLineColor(kRed);
                maxLineCD->Draw();

                canvas->cd(2);
                hadronMomentumFDHists.at(hadron)->Draw();

                TLine *minLineFD = new TLine(MIN_FD, 0, MIN_FD, hadronMomentumFDHists.at(hadron)->GetMaximum());
                minLineFD->SetLineColor(kRed);
                minLineFD->Draw();

                TLine *maxLineFD = new TLine(MAX, 0, MAX, hadronMomentumFDHists.at(hadron)->GetMaximum());
                maxLineFD->SetLineColor(kRed);
                maxLineFD->Draw();

                canvas->Update();
                canvas->SaveAs("CUT_Pi_Plus_Momentum.png");

                delete minLineCD;
                delete maxLineCD;
                delete minLineFD;
                delete maxLineFD;
                delete canvas;
            }
        }
    }
private:
    const double MIN_FD = 0.4;
    const double MIN_CD = 0.2;
    const double MAX = BEAM.P();

    std::map<Pids, TH1*> hadronMomentumCDHists;
    std::map<Pids, TH1*> hadronMomentumFDHists;
};

class HadronBetaCut : public Cut<const DataBanks&, int> {
public:
    HadronBetaCut() {
        for (Pids hadron : CHARGED_HADRONS) {
            std::string title = "Hadron CD #beta_{" + std::to_string(hadron) + "}; #beta";
            std::string name = "hadron_CD_beta_" + std::to_string(hadron);
            hadronBetaCDHists[hadron] = new TH1F(name.c_str(), title.c_str(), 200, 0.1, 1.2);

            title = "Hadron FD #beta_{" + std::to_string(hadron) + "}; #beta";
            name = "hadron_FD_beta_" + std::to_string(hadron);
            hadronBetaFDHists[hadron] = new TH1F(name.c_str(), title.c_str(), 200, 0.3, 1.2);
        }
    }

    ~HadronBetaCut() {
        for (auto it = hadronBetaCDHists.begin(); it != hadronBetaCDHists.end(); it++) {
            delete it->second;
        }

        for (auto it = hadronBetaFDHists.begin(); it != hadronBetaFDHists.end(); it++) {
            delete it->second;
        }
    }

    bool operator()(const DataBanks& banks, int i) override {
        const hipo::bank* PART = &banks.PART;

        int status = abs(PART->getInt("status", i));
        double beta = PART->getFloat("beta", i);
        Pids pid = Pids(PART->getInt("pid", i));

        if (status >= 2000 && status < 4000) {
            hadronBetaFDHists[pid]->Fill(beta);
            return beta > MIN_FD && beta < MAX_FD;
        }
        hadronBetaCDHists[pid]->Fill(beta);
        return beta > MIN_CD && beta < MAX_CD;
    }

    void check() const override {
        for (Pids hadron : CHARGED_HADRONS) {
            if (hadron == PI_PLUS) {
                TCanvas* canvas = new TCanvas("canvas_pi_plus_beta", "Histogrtams", 1200, 450);
                canvas->Divide(2, 1);

                canvas->cd(1);
                hadronBetaCDHists.at(hadron)->Draw();

                TLine *minLineCD = new TLine(MIN_CD, 0, MIN_CD, hadronBetaCDHists.at(hadron)->GetMaximum());
                minLineCD->SetLineColor(kRed);
                minLineCD->Draw();

                TLine *maxLineCD = new TLine(MAX_CD, 0, MAX_CD, hadronBetaCDHists.at(hadron)->GetMaximum());
                maxLineCD->SetLineColor(kRed);
                maxLineCD->Draw();

                canvas->cd(2);
                hadronBetaFDHists.at(hadron)->Draw();

                TLine *minLineFD = new TLine(MIN_FD, 0, MIN_FD, hadronBetaFDHists.at(hadron)->GetMaximum());
                minLineFD->SetLineColor(kRed);
                minLineFD->Draw();

                TLine *maxLineFD = new TLine(MAX_FD, 0, MAX_FD, hadronBetaFDHists.at(hadron)->GetMaximum());
                maxLineFD->SetLineColor(kRed);
                maxLineFD->Draw();

                canvas->Update();
                canvas->SaveAs("CUT_Pi_Plus_Beta.png");

                delete minLineCD;
                delete maxLineCD;
                delete minLineFD;
                delete maxLineFD;
                delete canvas;
            }
        }
    }
private:
    const double MIN_FD = 0.4;
    const double MAX_FD = 1.1;
    const double MIN_CD = 0.2;
    const double MAX_CD = 1.1;

    std::map<Pids, TH1*> hadronBetaCDHists;
    std::map<Pids, TH1*> hadronBetaFDHists;
};

class HadronZVertexCut : public Cut<const DataBanks&, int> {
public:
    HadronZVertexCut() {
        for (Pids hadron : CHARGED_HADRONS) {
            std::string title = "Hadron {" + std::to_string(hadron) + "} Z Vertex; v_{z}, cm";
            std::string name = "hadron_" + std::to_string(hadron) + "_z_vertex";
            hadronZVertexHists[hadron] = new TH1F(name.c_str(), title.c_str(), 200, -12.0, 4.0);
        }
    }

    ~HadronZVertexCut() {
        for (auto it = hadronZVertexHists.begin(); it != hadronZVertexHists.end(); it++) {
            delete it->second;
        }
    }

    bool operator()(const DataBanks& banks, int i) override {
        const hipo::bank* PART = &banks.PART;
    
        double vz = PART->getFloat("vz", i);

        hadronZVertexHists[Pids(PART->getInt("pid", i))]->Fill(vz);

        return vz > MIN && vz < MAX;
    }

    void check() const override {
        for (Pids hadron : CHARGED_HADRONS) {
            if (hadron == PI_PLUS) {
                TCanvas* canvas = new TCanvas("canvas_pi_plus_z_vertex", "Histogrtams", 600, 450);
                canvas->Divide(1, 1);

                canvas->cd(1);
                hadronZVertexHists.at(hadron)->Draw();

                TLine *minLine = new TLine(MIN, 0, MIN, hadronZVertexHists.at(hadron)->GetMaximum());
                minLine->SetLineColor(kRed);
                minLine->Draw();

                TLine *maxLine = new TLine(MAX, 0, MAX, hadronZVertexHists.at(hadron)->GetMaximum());
                maxLine->SetLineColor(kRed);
                maxLine->Draw();

                canvas->Update();
                canvas->SaveAs("CUT_Pi_Plus_Z_Vertex.png");

                delete minLine;
                delete maxLine;
                delete canvas;
            }
        }
    }
private:
    const double MIN = -10.0;
    const double MAX = 2.0;

    std::map<Pids, TH1*> hadronZVertexHists;
};

class HadronDCFiducialCut : public Cut<const DataBanks&, int> {
public:
    HadronDCFiducialCut() {
        for (Pids part : CHARGED_HADRONS_AND_ELECTRON) {
            for (auto& layer : LAYERS) {
                std::string title = "DC Particle " + std::to_string(part) + " y Vs x in L" + std::to_string(layer) + "; x, cm; y, cm";
                std::string name = "DC_part_" + std::to_string(part) + "_layer_" + std::to_string(layer);
                particleDCFiducialHists[part][layer] = new TH2F(name.c_str(), title.c_str(), 200, 20.0, 350.0, 200, -150.0, 150.0);
            }
        }
    }

    ~HadronDCFiducialCut() {
        for (auto it = particleDCFiducialHists.begin(); it != particleDCFiducialHists.end(); it++) {
            for (auto hist = it->second.begin(); hist != it->second.end(); hist++) {
                delete hist->second;
            }
        }
    }

    bool operator()(const DataBanks& banks, int i) override {
        const hipo::bank* TRAJ = &banks.TRAJ;
        const hipo::bank* TRACK = &banks.TRACK;
        const hipo::bank* PART = &banks.PART;

        const int DETECTOR = 6;

        for (int j = 0; j < TRAJ->getRows(); j++) {
            if (TRAJ->getInt("pindex", j) == i && TRAJ->getInt("detector", j) == DETECTOR) {
                for (int k = 0; k < TRACK->getRows(); k++) {
                    if (TRACK->getInt("pindex", k) == i && TRACK->getInt("detector", k) == DETECTOR) {
                        double* params = DCFiducialCutParams(TRAJ->getInt("layer", j), PART->getInt("charge", i));

                        double x0 = TRAJ->getFloat("x", j);
                        double y0 = TRAJ->getFloat("y", j);
                        int sector = findSector(x0, y0);

                        std::pair<double, double> coordinate = rotatePoint(x0, y0, sector);

                        double x = coordinate.first;
                        double y = coordinate.second;

                        particleDCFiducialHists[Pids(PART->getInt("pid", i))][TRAJ->getInt("layer", j)]->Fill(x, y);

                        if (
                            y > params[0] * x + params[1] ||
                            y < params[2] * x + params[3] ||
                            x < params[4]
                        ) {
                            delete params;
                            return false;
                        }

                        delete params;
                    }
                }
            }
        }
        return true;
    }

    void check() const override {
        for (Pids part : CHARGED_HADRONS_AND_ELECTRON) {
            if (part == PI_PLUS || part == ELECTRON) {
                std::string canvasName = "canvas_part_" + std::to_string(part) + "_drift_camera";
                TCanvas* canvas = new TCanvas(canvasName.c_str(), "Histogrtams", 1800, 450);
                canvas->Divide(3, 1);

                vector<TLine*> upLines, lowLines, edgeLines;
                upLines.resize(SECTORS);
                lowLines.resize(SECTORS);
                edgeLines.resize(SECTORS);

                for (int i = 0; i < SECTORS; i++) {
                    int layer = LAYERS[i];
                    canvas->cd(i + 1);
                    particleDCFiducialHists.at(part).at(layer)->Draw("colz");

                    double* params = DCFiducialCutParams(layer, partCharge(part));

                    upLines[i] = new TLine(params[4], params[0] * params[4] + params[1], 250, params[0] * 250 + params[1]);
                    lowLines[i] = new TLine(params[4], params[2] * params[4] + params[3], 250, params[2] * 250 + params[3]);
                    edgeLines[i] = new TLine(params[4], params[0] * params[4] + params[1], params[4], params[2] * params[4] + params[3]);

                    upLines[i]->SetLineColor(kRed);
                    lowLines[i]->SetLineColor(kRed);
                    edgeLines[i]->SetLineColor(kRed);

                    upLines[i]->Draw();
                    lowLines[i]->Draw();
                    edgeLines[i]->Draw();

                    delete params;
                }

                canvas->Update();

                std::string saveName = "CUT_Particle_" + std::to_string(part) + "_Drift_Camera.png";
                canvas->SaveAs(saveName.c_str());

                for (int i = 0; i < SECTORS; i++) {
                    delete upLines[i];
                    delete lowLines[i];
                    delete edgeLines[i];
                }

                delete canvas;
            }
        }
    }
private:
    static const int R1 = 6;
    static const int R2 = 18;
    static const int R3 = 36;
    static const int SECTORS = 3;
    const int LAYERS[SECTORS] = { R1, R2, R3 };

    std::map<Pids, std::map<int, TH2F*>> particleDCFiducialHists;

    int findSector(double x, double y) const {
        double angle = atan2(y, x);
        if (angle < 0) {
            angle += 2 * M_PI;
        }
        
        int sector = (int)((angle + M_PI / 6) / (M_PI / 3)) % 6 + 1;
        return sector;
    }

    double* DCFiducialCutParams(int layer, int charge) const {
        double* params = new double[5];

        int sgn = charge > 0 ? 1 : charge < 0 ? -1 : 0;

        switch (layer * sgn)
        {
        case R1:
            return new double[5] {0.610, -12.720, -0.604, 12.159, 38.02};
        case R2:
            return new double[5] {0.573, -13.949, -0.569, 13.891, 54.88};
        case R3:
            return new double[5] {0.527, -11.998, -0.530, 13.372, 49.0};
        case -R1:
            return new double[5] {0.556, -6.878, -0.560, 7.482, 24.052};
        case -R2:
            return new double[5] {0.578, -13.898, -0.577, 14.851, 39.705};
        case -R3:
            return new double[5] {0.591, -27.459, -0.588, 26.912, 77.755};
        default:
            return new double[0];
        }
    }
};

std::pair<double, double> rotatePoint(double x, double y, int sector) {
    const double ANGLE = M_PI / 3;
    double rotate = -ANGLE * (sector - 1);
    return std::make_pair(
        x * cos(rotate) - y * sin(rotate),
        x * sin(rotate) + y * cos(rotate)
    );
}

TGraph* getGraph(function<double(double)> func, double xmin, double xmax, int points) {
        TGraph* graph = new TGraph();
        double step = (xmax - xmin) / points;

        for (int i = 0; i < points; i++) {
            double x = xmin + i * step;
            graph->SetPoint(i, x, func(x));
        }

        graph->SetLineColor(kRed);
        return graph;
    }
