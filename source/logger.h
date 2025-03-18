#pragma once

#include <iostream>
#include <string>

#include "dataBanks.h"

using namespace std;

enum LoggerMode {
    QUIET = 0,
    WORK = 1,
    DEBUG = 2
};

class Logger {
public:
    Logger(LoggerMode mode=WORK) : numberEvent(0), mode(mode) {}

    unsigned int getNumberEvent() {
        return numberEvent;
    }

    LoggerMode getMode() {
        return mode;
    }

    void info(const DataBanks &banks) {
        const int REPEAT = 100000;
        numberEvent++;

        switch (mode)
        {
        case QUIET:
            break;
        case WORK:
            if (numberEvent % REPEAT == 0) {
                cout << "Passed Event: " << numberEvent << endl;
            }
            break;
        case DEBUG:
            if (numberEvent % REPEAT == 0) {
                cout << "Passed Event: " << numberEvent << endl;
                /**
                for (auto& bank : banks) {
                    bank->show();
                }
                */
            }
            break;
        default:
            break;
        }
    }

private:
    unsigned int numberEvent;
    LoggerMode mode;
};