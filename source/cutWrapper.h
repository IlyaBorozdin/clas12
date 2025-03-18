#pragma once

#include <functional>
#include <vector>

#include "dataBanks.h"
#include "constants.h"
#include "rg-k.h"

using CutWrapper = Cut<const DataBanks&>;

class CommonCuts : public CutWrapper {
public:
    using FunctionType = Cut<const DataBanks&>;

    ~CommonCuts() {
        for (auto elem : funcs) {
            delete elem;
        }
    }
    CommonCuts(std::vector<FunctionType*> funcs) : funcs(funcs) {}

    bool operator()(const DataBanks& banks) override {
        for (auto& cut : funcs) {
            if (!(*cut)(banks)) {
                return false;
            }
        }
        return true;
    }

    void check() const override {
        for (auto& cut : funcs) {
            cut->check();
        }
    }

protected:
    std::vector<FunctionType*> funcs;
};

class HadronCuts : public CutWrapper {
public:
    using FunctionType = Cut<const DataBanks&, int>;

    ~HadronCuts() {
        for (auto elem : funcs) {
            delete elem;
        }
    }
    HadronCuts(std::vector<FunctionType*> funcs) : funcs(funcs) {}

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        for (int i = 0; i < PART->getRows(); i++) {
            if (isChargedHadron(Pids(PART->getInt("pid", i)))) {
                for (auto& cut : funcs) {
                    if (!(*cut)(banks, i)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    void check() const override {
        for (auto& cut : funcs) {
            cut->check();
        }
    }

protected:
    std::vector<FunctionType*> funcs;
};

class HadronAndElectronCuts : public HadronCuts {
public:
    using FunctionType = Cut<const DataBanks&, int>;

    HadronAndElectronCuts(std::vector<FunctionType*> funcs) : HadronCuts(funcs) {}

    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;

        for (int i = 0; i < PART->getRows(); i++) {
            if (isChargedHadronOrElectron(Pids(PART->getInt("pid", i)))) {
                for (auto& cut : funcs) {
                    if (!(*cut)(banks, i)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
};
