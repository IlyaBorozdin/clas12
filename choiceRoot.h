#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "source/dataAnalysis.h"
#include "topologyCut.h"

using ElectronCuts = CommonCuts;

// ================================================================
// === БАЗОВЫЙ КЛАСС ==============================================
// ================================================================

class ChoiceRootBase : public HipoDataAnalysis {
public:
    ChoiceRootBase(const char* inputFileName,
                   const char* outFileName,
                   bool appendMode = false,
                   const char* treeName = "ExpData")
        : HipoDataAnalysis(readInputFile(inputFileName)),
          file(outFileName, appendMode ? "UPDATE" : "RECREATE"),
          tree(treeName, "Tree with physics data") {
        // ВНИМАНИЕ: больше не вызываем setDefaultCuts() здесь!
        // Это важно, потому что виртуальные методы в конструкторах не виртуальны.
        // Вызываем только setupBranches (безопасно).
        setupBranches();
    }

    virtual ~ChoiceRootBase() = default;

    void analysisEvent(const DataBanks& banks) override {
        fillFromBanks(banks);
        tree.Fill();
    }

    void results() override {
        file.cd();
        tree.Write();
        file.Close();
    }

protected:
    TFile file;
    TTree tree;

    double e_px{}, e_py{}, e_pz{};
    double pi_px{}, pi_py{}, pi_pz{};

    virtual void setupBranches() {
        tree.Branch("e_px", &e_px, "e_px/D");
        tree.Branch("e_py", &e_py, "e_py/D");
        tree.Branch("e_pz", &e_pz, "e_pz/D");
        tree.Branch("pi_px", &pi_px, "pi_px/D");
        tree.Branch("pi_py", &pi_py, "pi_py/D");
        tree.Branch("pi_pz", &pi_pz, "pi_pz/D");
    }

    virtual void fillFromBanks(const DataBanks& banks) = 0;
    virtual void setDefaultCuts() {}  // ничего по умолчанию

    static std::vector<std::string> readInputFile(const char* inputFileName) {
        std::vector<std::string> dataFileNames;
        std::ifstream inputFile(inputFileName);
        if (!inputFile.is_open()) {
            std::cerr << "Incorrect File: " << inputFileName << std::endl;
            return dataFileNames;
        }
        std::string line;
        while (std::getline(inputFile, line)) {
            line.erase(std::find_if(line.rbegin(), line.rend(),
                                    [](unsigned char ch) { return !std::isspace(ch); })
                           .base(),
                       line.end());
            if (!line.empty())
                dataFileNames.push_back(line);
        }
        return dataFileNames;
    }

    // --- общий метод, создающий стандартные наборы Cuts ---
    static std::vector<CutWrapper*> makeStandardCuts() {
        std::vector<Cut<const DataBanks&>*> electronCuts = {
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

        std::vector<Cut<const DataBanks&, int>*> hadronCuts = {
            new HadronMomentumCut(),
            new HadronBetaCut(),
            new HadronZVertexCut(),
        };

        std::vector<Cut<const DataBanks&, int>*> hadronAndElectronCuts = {
            new HadronDCFiducialCut(),
        };

        return {
            new ElectronCuts(electronCuts),
            new TopologyCut(),
            new HadronCuts(hadronCuts),
            new HadronAndElectronCuts(hadronAndElectronCuts),
        };
    }
};

// ================================================================
// === ChoiceRoot — ДАННЫЕ С РЕАЛЬНОГО ЭКСПЕРИМЕНТА ===============
// ================================================================

class ChoiceRoot : public ChoiceRootBase {
public:
    // явный конструктор, чтобы после полной инициализации вызвать setDefaultCuts()
    ChoiceRoot(const char* inputFileName,
               const char* outFileName,
               bool appendMode = false,
               const char* treeName = "ExpData")
        : ChoiceRootBase(inputFileName, outFileName, appendMode, treeName) {
        // теперь виртуальный setDefaultCuts() вызывает именно ChoiceRoot::setDefaultCuts()
        setDefaultCuts();
    }

protected:
    void setDefaultCuts() override {
        this->setStandartCuts(makeStandardCuts());
    }

    void fillFromBanks(const DataBanks& banks) override {
        const auto& PART = banks.PART;
        for (int i = 0; i < PART.getRows(); i++) {
            switch (PART.getInt("pid", i)) {
                case ELECTRON:
                    e_px = PART.getFloat("px", i);
                    e_py = PART.getFloat("py", i);
                    e_pz = PART.getFloat("pz", i);
                    break;
                case PI_PLUS:
                    pi_px = PART.getFloat("px", i);
                    pi_py = PART.getFloat("py", i);
                    pi_pz = PART.getFloat("pz", i);
                    break;
                default:
                    break;
            }
        }
    }
};

// ================================================================
// === ChoiceRootWeight — ВЕСОВЫЕ СОБЫТИЯ =========================
// ================================================================

class ChoiceRootWeight : public ChoiceRootBase {
public:
    using ChoiceRootBase::ChoiceRootBase;

    void analysisEvent(const DataBanks& banks) override {
        fillFromBanks(banks);
        correctWeight(weight);
        tree.Fill();
    }

protected:
    double weight{};

    void setupBranches() override {
        ChoiceRootBase::setupBranches();
        tree.Branch("weight", &weight, "weight/D");
    }

    virtual void correctWeight(double& weight) {}
};

// ================================================================
// === ChoiceRootSim — СИМУЛЯЦИЯ С ТЕМИ ЖЕ CUTS ===================
// ================================================================

class ChoiceRootSim : public ChoiceRootWeight {
public:
    // явный конструктор — вызываем setDefaultCuts() здесь
    ChoiceRootSim(const char* inputFileName,
                  const char* outFileName,
                  bool appendMode = false,
                  const char* treeName = "ExpData")
        : ChoiceRootWeight(inputFileName, outFileName, appendMode, treeName) {
        // ВАЖНО: этот вызов теперь сработает как полиморфный — будет вызван ChoiceRootSim::setDefaultCuts
        setDefaultCuts();
    }

protected:
    double trueN_px{}, trueN_py{}, trueN_pz{};
    double trueG_px{}, trueG_py{}, trueG_pz{};
    double trueE_px{}, trueE_py{}, trueE_pz{};
    double trueP_px{}, trueP_py{}, trueP_pz{};

    void setupBranches() override {
        ChoiceRootWeight::setupBranches();
        tree.Branch("trueE_px", &trueE_px, "trueE_px/D");
        tree.Branch("trueE_py", &trueE_py, "trueE_py/D");
        tree.Branch("trueE_pz", &trueE_pz, "trueE_pz/D");
        tree.Branch("trueP_px", &trueP_px, "trueP_px/D");
        tree.Branch("trueP_py", &trueP_py, "trueP_py/D");
        tree.Branch("trueP_pz", &trueP_pz, "trueP_pz/D");
        tree.Branch("trueN_px", &trueN_px, "trueN_px/D");
        tree.Branch("trueN_py", &trueN_py, "trueN_py/D");
        tree.Branch("trueN_pz", &trueN_pz, "trueN_pz/D");
        tree.Branch("trueG_px", &trueG_px, "trueG_px/D");
        tree.Branch("trueG_py", &trueG_py, "trueG_py/D");
        tree.Branch("trueG_pz", &trueG_pz, "trueG_pz/D");
    }

    void setDefaultCuts() override {
        // выбор тех же стандартных наборов, что в ChoiceRoot
        this->setStandartCuts(makeStandardCuts());
    }

    void fillFromBanks(const DataBanks& banks) override {
        const auto& PART = banks.PART;
        const auto& LUND = banks.LUND;
        const auto& EVENT = banks.EVENT;

        for (int i = 0; i < PART.getRows(); i++) {
            switch (PART.getInt("pid", i)) {
                case ELECTRON:
                    e_px = PART.getFloat("px", i);
                    e_py = PART.getFloat("py", i);
                    e_pz = PART.getFloat("pz", i);
                    break;
                case PI_PLUS:
                    pi_px = PART.getFloat("px", i);
                    pi_py = PART.getFloat("py", i);
                    pi_pz = PART.getFloat("pz", i);
                    break;
                default:
                    break;
            }
        }

        for (int i = 0; i < LUND.getRows(); i++) {
            switch (LUND.getInt("pid", i)) {
                case ELECTRON:
                    trueE_px = LUND.getFloat("px", i);
                    trueE_py = LUND.getFloat("py", i);
                    trueE_pz = LUND.getFloat("pz", i);
                    break;
                case PI_PLUS:
                    trueP_px = LUND.getFloat("px", i);
                    trueP_py = LUND.getFloat("py", i);
                    trueP_pz = LUND.getFloat("pz", i);
                    break;
                case NEUTRON:
                    trueN_px = LUND.getFloat("px", i);
                    trueN_py = LUND.getFloat("py", i);
                    trueN_pz = LUND.getFloat("pz", i);
                    break;
                case PHOTON:
                    trueG_px = LUND.getFloat("px", i);
                    trueG_py = LUND.getFloat("py", i);
                    trueG_pz = LUND.getFloat("pz", i);
                    break;
                default:
                    break;
            }
        }

        weight = EVENT.getFloat("weight", 0);
    }

    void correctWeight(double& weight) override {
        TLorentzVector electron, piPlus;
        electron.SetXYZM(trueE_px, trueE_py, trueE_pz, ELECTRON_MASS);
        piPlus.SetXYZM(trueP_px, trueP_py, trueP_pz, CHARGED_PI_MASS);

        TLorentzVector gamma = BEAM - electron;

        TLorentzVector piPlusCM = piPlus;
        TVector3 betaLab(0.0, 0.0, -gamma.P() / (gamma.E() + PROTON_MASS));
        piPlusCM.RotateZ(-electron.Phi());
        piPlusCM.RotateY(gamma.Theta());
        piPlusCM.Boost(betaLab);
        piPlusCM.SetXYZM(-piPlusCM.Px(), -piPlusCM.Py(), -piPlusCM.Pz(), CHARGED_PI_MASS);

        weight *= TMath::Sin(piPlusCM.Theta());
    }
};

// ================================================================
// === ChoiceRootLund — РАБОТА С LUND-БАНКОМ =======================
// ================================================================

class ChoiceRootLund : public ChoiceRootWeight {
public:
    using ChoiceRootWeight::ChoiceRootWeight;

protected:
    double n_px{}, n_py{}, n_pz{};
    double g_px{}, g_py{}, g_pz{};

    void setupBranches() override {
        ChoiceRootWeight::setupBranches();
        tree.Branch("n_px", &n_px, "n_px/D");
        tree.Branch("n_py", &n_py, "n_py/D");
        tree.Branch("n_pz", &n_pz, "n_pz/D");
        tree.Branch("g_px", &g_px, "g_px/D");
        tree.Branch("g_py", &g_py, "g_py/D");
        tree.Branch("g_pz", &g_pz, "g_pz/D");
    }

    void fillFromBanks(const DataBanks& banks) override {
        const auto& LUND = banks.LUND;
        const auto& EVENT = banks.EVENT;

        for (int i = 0; i < LUND.getRows(); i++) {
            switch (LUND.getInt("pid", i)) {
                case ELECTRON:
                    e_px = LUND.getFloat("px", i);
                    e_py = LUND.getFloat("py", i);
                    e_pz = LUND.getFloat("pz", i);
                    break;
                case PI_PLUS:
                    pi_px = LUND.getFloat("px", i);
                    pi_py = LUND.getFloat("py", i);
                    pi_pz = LUND.getFloat("pz", i);
                    break;
                case NEUTRON:
                    n_px = LUND.getFloat("px", i);
                    n_py = LUND.getFloat("py", i);
                    n_pz = LUND.getFloat("pz", i);
                    break;
                case PHOTON:
                    g_px = LUND.getFloat("px", i);
                    g_py = LUND.getFloat("py", i);
                    g_pz = LUND.getFloat("pz", i);
                    break;
                default:
                    break;
            }
        }

        weight = EVENT.getFloat("weight", 0);
    }

    void correctWeight(double& weight) override {
        TLorentzVector electron, piPlus;
        electron.SetXYZM(e_px, e_py, e_pz, ELECTRON_MASS);
        piPlus.SetXYZM(pi_px, pi_py, pi_pz, CHARGED_PI_MASS);

        TLorentzVector gamma = BEAM - electron;

        TLorentzVector piPlusCM = piPlus;
        TVector3 betaLab(0.0, 0.0, -gamma.P() / (gamma.E() + PROTON_MASS));
        piPlusCM.RotateZ(-electron.Phi());
        piPlusCM.RotateY(gamma.Theta());
        piPlusCM.Boost(betaLab);
        piPlusCM.SetXYZM(-piPlusCM.Px(), -piPlusCM.Py(), -piPlusCM.Pz(), CHARGED_PI_MASS);

        weight *= TMath::Sin(piPlusCM.Theta());
    }
};
