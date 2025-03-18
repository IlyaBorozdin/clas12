#pragma once

#include <iterator>

#include "clas12reader.h"

using namespace clas12;

class DataBanks {
public:
    hipo::bank PART;
    hipo::bank SCINT;
    hipo::bank ECAL;
    hipo::bank TRAJ;
    hipo::bank TRACK;

    auto begin() const { return std::begin(banks); }
    auto end() const { return std::end(banks); }

    void getSchema(hipo::dictionary &factory) {
        PART = factory.hasSchema("REC::Particle") ? hipo::bank(factory.getSchema("REC::Particle")) : hipo::bank();
        SCINT = factory.hasSchema("REC::Scintillator") ? hipo::bank(factory.getSchema("REC::Scintillator")) : hipo::bank();
        ECAL = factory.hasSchema("REC::Calorimeter") ? hipo::bank(factory.getSchema("REC::Calorimeter")) : hipo::bank();
        TRAJ = factory.hasSchema("REC::Traj") ? hipo::bank(factory.getSchema("REC::Traj")) : hipo::bank();
        TRACK = factory.hasSchema("REC::Track") ? hipo::bank(factory.getSchema("REC::Track")) : hipo::bank();
    }

    void getStructure(hipo::event &event) {
        for (auto& bank : *this) {
            event.getStructure(*bank);
        }
    }

private:
    hipo::bank* banks[5] = {&PART, &SCINT, &ECAL, &TRAJ, &TRACK};
};
