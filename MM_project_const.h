#pragma once

const double W_MIN = 1.1;
const double W_MAX = 2.2;

const int NUMBER_Q2 = 7;
const int NUMBER_W = 44;
const int NUMBER_COS_THETA = 10;
const int NUMBER_PHI = 10;

const double STEP_W = (W_MAX - W_MIN) / NUMBER_W;
const double STEPS_Q2[14] = { 0.4, 0.6, 0.6, 0.8, 0.8, 1.0, 1.0, 1.4, 1.4, 1.8, 1.8, 2.6, 2.6, 3.5 };
const double STEP_COS_THETA = 2.0 / NUMBER_COS_THETA;
const double STEP_PHI = 2 * M_PI / NUMBER_PHI;

void convertNumbersToVariable(int& i, int& j, int& k, int& l, double& q2, double& w, double& cos_theta, double& phi) {
        q2 = (STEPS_Q2[2 * i] + STEPS_Q2[2 * i + 1]) / 2;
        w = (0.5 + j) * STEP_W + W_MIN;
        cos_theta = (0.5 + k) * STEP_COS_THETA - 1;
        phi = (0.5 + l) * STEP_PHI;
    }

const int OUTPUT_ENTRIES = 79948394;
const int OUTPUTSIM_ENTRIES = 62842196; // 125921656;
const int OUTPUTLUND_ENTRIES = 249265000; // 125920346;

const char* GEN_BASE_PATH = std::getenv("CLAS12_GEN_BASE");
const char* HEAVY_BASE_PATH = std::getenv("CLAS12_HEAVY_BASE");
