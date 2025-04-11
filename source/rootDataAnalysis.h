#pragma once

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

class RootDataAnalysis {
public:
    RootDataAnalysis(const std::string& rootFileName) 
        : file(rootFileName.c_str(), "READ"), tree(nullptr), numberEvent(0) {
        
        if (!file.IsOpen()) {
            std::cerr << "Error opening ROOT file: " << rootFileName << std::endl;
            return;
        }

        file.GetObject("ExpData", tree);
        if (!tree) {
            std::cerr << "Error: Tree 'ExpData' not found in file " << rootFileName << std::endl;
            return;
        }

        tree->SetBranchAddress("e_px", &e_px);
        tree->SetBranchAddress("e_py", &e_py);
        tree->SetBranchAddress("e_pz", &e_pz);
        tree->SetBranchAddress("pi_px", &pi_px);
        tree->SetBranchAddress("pi_py", &pi_py);
        tree->SetBranchAddress("pi_pz", &pi_pz);
    }

    virtual ~RootDataAnalysis() {
        if (file.IsOpen()) {
            file.Close();
        }
    }

    virtual void analysisEvent() = 0;
    virtual void results() = 0;

    void analysisCycle() {
        if (!tree) return;

        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            tree->GetEntry(i);
            numberEvent++;

            if (numberEvent % REPEAT == 0) {
                std::cout << "Processed Event: " << numberEvent << std::endl;
            }

            analysisEvent();
        }

        results();
    }

protected:
    TFile file;
    TTree* tree;
    unsigned int numberEvent;
    const int REPEAT = 1000000;

    double e_px, e_py, e_pz;
    double pi_px, pi_py, pi_pz;
};
