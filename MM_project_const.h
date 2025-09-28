#pragma once

const double W_MIN = 1.1;
// const double W_MAX = 2.0;
const double W_MAX = 2.2;

// const int NUMBER_Q2 = 4;
const int NUMBER_Q2 = 7;
// const int NUMBER_W = 18;
// const int NUMBER_W = 36;
const int NUMBER_W = 22;
const int NUMBER_COS_THETA = 10;
const int NUMBER_PHI = 10;

const double STEP_W = (W_MAX - W_MIN) / NUMBER_W;
// const double STEPS_Q2[8] = { 0.4, 0.6, 0.6, 1.0, 1.0, 2.0, 2.0, 3.5 };
const double STEPS_Q2[14] = { 0.4, 0.6, 0.6, 0.8, 0.8, 1.0, 1.0, 1.4, 1.4, 1.8, 1.8, 2.6, 2.6, 3.5 };
const double STEP_COS_THETA = 2.0 / NUMBER_COS_THETA;
const double STEP_PHI = 2 * M_PI / NUMBER_PHI;

const int OUTPUT_ENTRIES = 79948394;
const int OUTPUTSIM_ENTRIES = 125921656;
const int OUTPUTLUND_ENTRIES = 125920346;