#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <iomanip>

#include "source/dataAnalysis.h"

/**
 * @brief Промежуточный класс-наследник HipoDataAnalysis.
 * Отвечает за управление выходным файлом отчета и хранение катов.
 */
class HipoStatsBase : public HipoDataAnalysis {
public:
    /**
     * @param hipoList Путь к .txt со списком .hipo файлов
     * @param outputFileName Путь к файлу, куда запишем статистику
     */
    HipoStatsBase(const char* hipoList, const char* outputFileName)
        : HipoDataAnalysis(hipoList), outputName(outputFileName) {}

    virtual ~HipoStatsBase() = default;

    /**
     * @brief Позволяет установить каты снаружи (через обертку Step)
     */
    void setCuts(const std::vector<CutWrapper*>& cuts) {
        this->setStandartCuts(cuts);
    }

    /**
     * @brief Реализация обязательного метода HipoDataAnalysis.
     * Вызывается автоматически в конце analysisCycle().
     */
    void results() override {
        std::ofstream outFile(outputName);
        if (!outFile.is_open()) {
            std::cerr << "[HipoStatsBase] ERROR: Could not open output file: " << outputName << std::endl;
            return;
        }

        writeHeader(outFile);
        writeReport(outFile);
        writeFooter(outFile);

        outFile.close();
    }

protected:
    std::string outputName;

    /**
     * @brief Чисто виртуальный метод для реализации в конкретных статистиках.
     * Именно здесь ты будешь писать os << "My Stat: " << value;
     */
    virtual void writeReport(std::ostream& os) = 0;

private:
    void writeHeader(std::ostream& os) {
        os << "====================================================\n";
        os << "       HIPO DATA ANALYSIS STATISTICS REPORT        \n";
        os << "====================================================\n\n";
    }

    void writeFooter(std::ostream& os) {
        os << "\n====================================================\n";
        os << "End of Report.\n";
    }
};

/**
 * @brief Инструмент для базового мониторинга состава частиц в HIPO файлах.
 */
class SimpleParticleStats : public HipoStatsBase {
public:
    // Используем конструктор базового класса
    using HipoStatsBase::HipoStatsBase;

    /**
     * @brief Анализ каждого отдельного события с проверкой на множественность.
     */
    void analysisEvent(const DataBanks& banks) override {
        totalEvents++;

        int nElectrons = 0;
        int nPiPlus = 0;

        const auto& PART = banks.PART;
        int rows = PART.getRows();

        // 1. Считаем количество интересующих нас частиц в конкретном событии
        for (int i = 0; i < rows; i++) {
            int pid = PART.getInt("pid", i);

            if (pid == ELECTRON) {
                nElectrons++;
            } else if (pid == PI_PLUS) {
                nPiPlus++;
            }
        }

        // 2. Условие: ровно один электрон
        if (nElectrons == 1) {
            electronEvents++; // Здесь храним события с n_e == 1
            
            // 3. Условие: ровно один электрон И ровно один пион
            if (nPiPlus == 1) {
                electronPiPlusEvents++; // Здесь n_e == 1 && n_pi == 1
            }
        }
    }

protected:
    /**
     * @brief Формирование итогового отчета.
     */
    void writeReport(std::ostream& os) override {
        // Рассчитываем доли (с проверкой на 0, чтобы не делить на ноль)
        double eFrac = (totalEvents > 0) ? (100.0 * electronEvents / totalEvents) : 0.0;
        double ePiFrac = (totalEvents > 0) ? (100.0 * electronPiPlusEvents / totalEvents) : 0.0;

        os << std::left << std::setw(40) << "1. Total events processed:" 
           << totalEvents << "\n";

        os << std::left << std::setw(40) << "2. Events with electron (e-):" 
           << electronEvents << " (" << std::fixed << std::setprecision(2) << eFrac << "%)\n";

        os << std::left << std::setw(40) << "3. Events with e- and pi+:" 
           << electronPiPlusEvents << " (" << std::fixed << std::setprecision(2) << ePiFrac << "%)\n";
           
        os << "----------------------------------------------------\n";
        if (totalEvents > 0 && electronEvents == 0) {
            os << "[!] WARNING: No electrons found in the data!\n";
        }
    }

private:
    // Счетчики
    unsigned long long totalEvents = 0;
    unsigned long long electronEvents = 0;
    unsigned long long electronPiPlusEvents = 0;
};