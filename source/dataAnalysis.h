#pragma once

#include <iostream>
// #include <vector>
// #include <functional>

// #include "clas12reader.h"

// #include "dataBanks.h"
#include "cutWrapper.h"
// #include "logger.h"

using namespace clas12;
using namespace std;

class HipoDataAnalysis {
public:
    HipoDataAnalysis(const vector<string> &dataFileNames) : dataFileNames(dataFileNames), numberEvent(0) {}

    virtual ~HipoDataAnalysis() {}

    virtual void analysisEvent(const DataBanks& banks) = 0;
    virtual void results() = 0;

    void analysisCycle() {
        for (const string& dataFileName : dataFileNames) {

            hipo::reader reader;
            reader.open(dataFileName.c_str());

            hipo::dictionary factory;
            reader.readDictionary(factory);

            DataBanks banks;
            hipo::event event;

            banks.getSchema(factory);

            while (reader.next()) {

                reader.read(event);
                banks.getStructure(event);

                numberEvent++;
                if (numberEvent % REPEAT == 0) {
                    std::cout << "Passed Event: " << numberEvent << std::endl;
                }

                try {
                    if (cutsPassed(banks)) {
                        analysisEvent(banks);
                    }
                }
                catch(const exception& e) {
                    std::cerr << e.what() << std::endl;
                }
            }
        }

        results();
    }

protected:
    void setStandartCuts(const vector<CutWrapper*> cuts) {
        standartCuts = cuts;
    }
    

private:
    const vector<string> dataFileNames;
    vector<CutWrapper*> standartCuts;
    unsigned int numberEvent;
    const int REPEAT = 100000;

    bool cutsPassed(const DataBanks& banks) const {
        for (auto& cut : standartCuts) {
            if (!(*cut)(banks)) {
                return false;
            }
        }
        return true;
    }
};
