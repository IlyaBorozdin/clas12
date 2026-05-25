#pragma once

#include <TLorentzVector.h>

enum Pids {
    ELECTRON = 11,
    POSITRON = -11,
    MU_MINUS = 13,
    MU_PLUS = -12,
    PROTON = 2212,
    NEUTRON = 2112,
    ANTIPROTON = -2212,
    PI0 = 111,
    PI_PLUS = 211,
    PI_MINUS = -211,
    K_PLUS = 321,
    K_MINUS = -321,
    DEUTRON = 45,
    PHOTON = 22,
    UNKNOWN = 0
};

const double ELECTRON_MASS = 0.0005;
const double MU_MASS = 0.1057;
const double PROTON_MASS = 0.9383;
const double NEUTRON_MASS = 0.9396;
const double NEUTRAL_PI_MASS = 0.1350;
const double CHARGED_PI_MASS = 0.1396;
const double CHARGED_K_MASS = 0.4937;
const double DEUTRON_MASS = 1.8756;
const double PHOTON_MASS = .0;

const double BEAM_ENERGY = 6.535;

const std::vector<Pids> CHARGED_HADRONS = { PROTON, ANTIPROTON, PI_PLUS, PI_MINUS, K_PLUS, K_MINUS };
const std::vector<Pids> CHARGED_HADRONS_AND_ELECTRON = { PROTON, ANTIPROTON, PI_PLUS, PI_MINUS, K_PLUS, K_MINUS, ELECTRON };
const std::vector<Pids> CHARGED_PARTICLES = { PROTON, ANTIPROTON, PI_PLUS, PI_MINUS, K_PLUS, K_MINUS, ELECTRON, POSITRON, MU_MINUS, MU_PLUS, DEUTRON };
const std::vector<Pids> CHARGED_PRODUCTS = { PROTON, ANTIPROTON, PI_PLUS, PI_MINUS, K_PLUS, K_MINUS, POSITRON, MU_MINUS, MU_PLUS, DEUTRON };
const std::vector<Pids> CHARGED_LEPTONS = { ELECTRON, POSITRON, MU_MINUS, MU_PLUS };

const TLorentzVector BEAM(.0, .0, BEAM_ENERGY, BEAM_ENERGY);
const TLorentzVector PROTON_TARGET(.0, .0, .0, PROTON_MASS);

double partMass(Pids pid) {
    switch (abs(pid))
    {
    case ELECTRON:
        return ELECTRON_MASS;
    case MU_MINUS:
        return MU_MASS;   
    case PROTON:
        return PROTON_MASS;   
    case NEUTRON:
        return NEUTRON_MASS;   
    case PI0:
        return NEUTRAL_PI_MASS;   
    case PI_PLUS:
        return CHARGED_PI_MASS;   
    case K_PLUS:
        return CHARGED_K_MASS;   
    case DEUTRON:
        return DEUTRON_MASS;
    case PHOTON:
        return PHOTON_MASS;
    default:
        return .0;
    }
}

int partCharge(Pids pid) {
    switch (pid)
    {
    case ELECTRON:
    case MU_MINUS:
    case ANTIPROTON:
    case PI_MINUS:
    case K_MINUS:
        return -1;
    case POSITRON:
    case MU_PLUS:
    case PROTON:
    case PI_PLUS:
    case K_PLUS:
    case DEUTRON:
        return 1;
    default:
        return 0;
    }
}

template <typename T>
bool isCharged(const std::vector<Pids>& chargedArray, T pid) {
    return std::find(chargedArray.begin(), chargedArray.end(), pid) != chargedArray.end();
}

bool isChargedHadron(Pids pid) {
    return isCharged(CHARGED_HADRONS, pid);
}

bool isChargedHadronOrElectron(Pids pid) {
    return isCharged(CHARGED_HADRONS_AND_ELECTRON, pid);
}

bool isChargedParticle(Pids pid) {
    return isCharged(CHARGED_PARTICLES, pid);
}

bool isChargedProduct(Pids pid) {
    return isCharged(CHARGED_PRODUCTS, pid);
}

bool isChargedLepton(Pids pid) {
    return isCharged(CHARGED_LEPTONS, pid);
}
