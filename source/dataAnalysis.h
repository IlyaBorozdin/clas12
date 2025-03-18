#pragma once

#include <iostream>
#include <vector>
#include <functional>

#include "clas12reader.h"

#include "dataBanks.h"
#include "logger.h"
#include "cutWrapper.h"

using namespace clas12;

class HipoDataAnalysis {
public:
    HipoDataAnalysis(const vector<string> &dataFileNames) : dataFileNames(dataFileNames), logger(WORK) {}

    virtual ~HipoDataAnalysis() {
        /*
        for (auto& cut : standartCuts) {
            delete cut;
        }
        standartCuts.clear();
        */
    }

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

                logger.info(banks);

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
        
        if (logger.getMode() == DEBUG) {
            for (auto& cut : standartCuts) {
                cut->check();
            }
        }
    }

protected:
    Logger logger;

    void setStandartCuts(const vector<CutWrapper*> cuts) {
        standartCuts = cuts;
    }
    

private:
    const vector<string> dataFileNames;
    vector<CutWrapper*> standartCuts;

    bool cutsPassed(const DataBanks& banks) const {
        for (auto& cut : standartCuts) {
            if (!(*cut)(banks)) {
                return false;
            }
        }
        return true;
    }
};
