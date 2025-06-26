#pragma once

#include <string>

const double W_MIN = 1.1;
const double W_MAX = 2.0;

const int NUMBER_W = 18;
// const int NUMBER_W = 36;
const int NUMBER_Q2 = 4;
const int NUMBER_COS_THETA = 10;
const int NUMBER_PHI = 10;

const double STEP_W = (W_MAX - W_MIN) / NUMBER_W;
const double STEPS_Q2[8] = { 0.4, 0.6, 0.6, 1.0, 1.0, 2.0, 2.0, 3.5 };
const double STEP_COS_THETA = 2.0 / NUMBER_COS_THETA;
const double STEP_PHI = 2 * M_PI / NUMBER_PHI;

const double MM_MIN = 0.8;
const double MM_MAX = 1.6;
// const int NUMBER_MM = 80;
const int NUMBER_MM = 50;
const double MM_BIN_WIDTH = (MM_MAX - MM_MIN) / NUMBER_MM;

const int ENTRIES_MIN = 600;

// std::string MM_FILE = "missingMasses_50.root";
// std::string MM_FILE = "missingMasses_25.root";
// std::string MM_FILE = "missingMasses_25_4D.root";
std::string MM_FILE = "missingMasses_50_4D.root";
std::string IMG_BUFF = "img/buff/";

std::string nameHistMM(int i, int j) {
    return "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1);
}

std::string nameHist4D_MM(int i, int j, int k, int l) {
    return "CELL " + std::to_string(i + 1) + "-" + std::to_string(j + 1) + "-" + std::to_string(k + 1) + "-" + std::to_string(l + 1);
}