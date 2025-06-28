#pragma once

#include "source/rg-k.h"

class TopologyCut : public Cut<const DataBanks&> {
public:
    bool operator()(const DataBanks& banks) override {
        const hipo::bank* PART = &banks.PART;
        int countPiPlus = 0;
        int countCharge = 0;

        for (int i = 0; i < PART->getRows(); i++) {
            int pid = PART->getInt("pid", i);

            countCharge += PART->getInt("charge", i) == 0 ? 0 : 1;

            if (pid == PI_PLUS) {
                countPiPlus++;
            }
        }

        if (
            countPiPlus != 1 ||
            countCharge != 2
        ) {
            return false;
        }

        return true;
    }

    void check() const override {}
};