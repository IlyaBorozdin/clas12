#pragma once

#include <fstream>
// #include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"

// #include <cmath>

#include "source/dataAnalysis.h"
// #include "source/cutWrapper.h"
// #include "source/constants.h"
// #include "source/rg-k.h"
#include "topologyCut.h"

using ElectronCuts = CommonCuts;

vector<string> readInputFile(const char* inputFileName);

class ChoiceRoot : public HipoDataAnalysis {
public:
    ChoiceRoot(const char* inputFileName,
               const char* outFileName,
               bool appendMode = false,
               const char* treeName = "ExpData")
        : HipoDataAnalysis(readInputFile(inputFileName)),
          file(outFileName, appendMode ? "UPDATE" : "RECREATE"),
          tree(treeName, "Tree with electron and pion data") {

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

        tree.Branch("e_px", &e_px, "e_px/D");
        tree.Branch("e_py", &e_py, "e_py/D");
        tree.Branch("e_pz", &e_pz, "e_pz/D");
        tree.Branch("pi_px", &pi_px, "pi_px/D");
        tree.Branch("pi_py", &pi_py, "pi_py/D");
        tree.Branch("pi_pz", &pi_pz, "pi_pz/D");
    }

    ~ChoiceRoot() {}

    void analysisEvent(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        for (int i = 0; i < PART->getRows(); i++) {
            switch (PART->getInt("pid", i))
            {
            case ELECTRON:
                e_px = PART->getFloat("px", i);
                e_py = PART->getFloat("py", i);
                e_pz = PART->getFloat("pz", i);
                break;
            case PI_PLUS:
                pi_px = PART->getFloat("px", i);
                pi_py = PART->getFloat("py", i);
                pi_pz = PART->getFloat("pz", i);
                break;
            default:
                break;
            }
        }

        tree.Fill();
    }

    void results() override {
        file.cd();
        tree.Write();
        file.Close();
    }

private:
    TFile file;
    TTree tree;

    double e_px, e_py, e_pz, pi_px, pi_py, pi_pz;
};

vector<string> readInputFile(const char* inputFileName) {
    vector<string> dataFileNames;

    std::ifstream inputFile(inputFileName);

        if (!inputFile.is_open()) {
            std::cerr << "Incorrect File: " << inputFileName << std::endl;
            return dataFileNames;
        }
        std::string line;
        while (std::getline(inputFile, line)) {
            line.erase(std::find_if(line.rbegin(), line.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), line.end());

            if (!line.empty()) {
                dataFileNames.push_back(line);
            }
        }

        inputFile.close();
        return dataFileNames;
}