#pragma once

#include <iostream>
#include <filesystem>
#include <csignal>
#include <csetjmp>
#include <chrono>
#include <iomanip>
#include <sstream>

#include "cutWrapper.h"

using namespace clas12;
using namespace std;

// ================================================================
// === Глобальный jump buffer для обработки SIGSEGV и т.п. ========
// ================================================================
static jmp_buf jump_buffer;
static string currentFileName = "undefined";
static int currentEventIndex = -1;

void signalHandler(int signal) {
    cerr << "\n[FATAL] Signal " << signal
         << " caught while processing file: " << currentFileName
         << ", event index: " << currentEventIndex
         << "\nContinuing with next file..." << endl;
    longjmp(jump_buffer, 1);
}

// ================================================================
// === HipoDataAnalysis ===========================================
// ================================================================

class HipoDataAnalysis {
public:
    explicit HipoDataAnalysis(const vector<string>& dataFileNames)
        : dataFileNames(dataFileNames), numberEvent(0) {
        installSignalHandlers();
    }

    virtual ~HipoDataAnalysis() = default;

    virtual void analysisEvent(const DataBanks& banks) = 0;
    virtual void results() = 0;

    // Главный цикл обработки
    void analysisCycle() {
        vector<string> failedFiles;
        auto startAll = chrono::steady_clock::now();

        if (standartCuts.empty()) {
            std::cerr << "[Warning] No cuts have been set! All events will be processed." << std::endl;
        }

        for (const string& dataFileName : dataFileNames) {
            currentFileName = dataFileName;
            currentEventIndex = -1;

            if (!filesystem::exists(dataFileName)) {
                log("File not found: " + dataFileName, "ERROR");
                failedFiles.push_back(dataFileName);
                continue;
            }

            log("Opening file: " + dataFileName, "INFO");

            // Установка точки возврата при краше
            if (setjmp(jump_buffer) != 0) {
                // После longjmp попадаем сюда
                failedFiles.push_back(dataFileName);
                continue;  // переходим к следующему файлу
            }

            try {
                processFile(dataFileName);
            } catch (const std::exception& e) {
                log("Exception while processing file " + dataFileName + ": " + e.what(), "ERROR");
                failedFiles.push_back(dataFileName);
            }
        }

        results();

        auto endAll = chrono::steady_clock::now();
        auto totalSec = chrono::duration_cast<chrono::seconds>(endAll - startAll).count();

        log("Analysis finished in " + to_string(totalSec) + " s", "INFO");

        if (!failedFiles.empty()) {
            cerr << "\n=== FAILED FILES ===" << endl;
            for (const auto& f : failedFiles)
                cerr << " - " << f << endl;
            cerr << "====================" << endl;
        }
    }

protected:
    void setStandartCuts(const vector<CutWrapper*>& cuts) {
        standartCuts = cuts;
    }

private:
    const vector<string> dataFileNames;
    vector<CutWrapper*> standartCuts;
    unsigned long long numberEvent = 0;
    const int REPEAT = 100000;

    // ==================== Основная логика ====================

    void processFile(const string& dataFileName) {
        hipo::reader reader;
        reader.open(dataFileName.c_str());

        hipo::dictionary factory;
        reader.readDictionary(factory);

        DataBanks banks;
        hipo::event event;
        banks.getSchema(factory);

        unsigned long long eventCount = 0;
        unsigned long long passedCount = 0; // <--- добавили счётчик прошедших cuts

        auto start = chrono::steady_clock::now();

        while (reader.next()) {
            currentEventIndex = eventCount;

            reader.read(event);
            banks.getStructure(event);
            ++numberEvent;
            ++eventCount;

            if (numberEvent % REPEAT == 0)
                log("Processed events: " + to_string(numberEvent), "INFO");

            try {
                if (cutsPassed(banks)) {
                    ++passedCount;             // <--- увеличиваем счётчик
                    analysisEvent(banks);      // если прошли все cuts — обрабатываем
                }
            } catch (const std::exception& e) {
                log("Error inside event " + to_string(eventCount) + ": " + e.what(), "WARN");
            }
        }

        auto end = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();

        stringstream ss;
        ss << "Finished file: " << dataFileName
        << " | passed: " << passedCount << " / " << eventCount  // <--- нужный формат
        << " | time: " << fixed << setprecision(2) << elapsed / 1000.0 << " s";

        log(ss.str(), "INFO");
    }


    bool cutsPassed(const DataBanks& banks) const {
        for (auto& cut : standartCuts) {
            if (!(*cut)(banks)) return false;
        }
        return true;
    }

    // ==================== Вспомогательные методы ====================

    static void installSignalHandlers() {
        signal(SIGSEGV, signalHandler);
        signal(SIGBUS, signalHandler);
        signal(SIGABRT, signalHandler);
    }

    static void log(const string& msg, const string& level = "INFO") {
        auto now = chrono::system_clock::now();
        time_t now_time = chrono::system_clock::to_time_t(now);
        tm* ptm = localtime(&now_time);
        cout << "[" << put_time(ptm, "%H:%M:%S") << "] "
             << "[" << setw(5) << left << level << "] " << msg << endl;
    }
};
