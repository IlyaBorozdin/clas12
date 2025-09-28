#pragma once

#include <TLorentzVector.h>
#include <TVector3.h>
#include <cmath>

#include "rootDataAnalysisStep.h"
#include "source/constants.h"
#include "MM_project_const.h"
#include "MM_project_utils.h"

/**
 * @brief Шаг анализа, сохраняющий события в четырёхмерный массив TTree'ов
 *        по бинам в переменных Q², W, cos(θ), φ.
 *
 * Наследует RootDataAnalysisStep и использует данные событий (импульсы e⁻ и π⁺)
 * для расчёта недостающей массы и углов в КС-системе, распределяя события
 * по соответствующим бинам. Каждое событие сохраняется в свой TTree.
 */
class BinnedTreeStep : public RootDataAnalysisStep {
public:
    /**
     * @param stepName       Название шага (используется в логах)
     * @param inputFileName  Входной ROOT-файл с деревом ExpData
     * @param outputFileName Имя выходного файла, куда будут записаны деревья
     */
    BinnedTreeStep(const std::string& stepName,
                   const std::string& inputFileName,
                   const std::string& outputFileName)
        : RootDataAnalysisStep(stepName, inputFileName, outputFileName) {}

    /**
     * @brief Удаляет все созданные TTree из памяти
     */
    /*
    ~BinnedTreeStep() override {
        for (auto& q2vec : trees)
            for (auto& wvec : q2vec)
                for (auto& ctvec : wvec)
                    for (auto* tree : ctvec)
                        delete tree;
    }
    */
    ~BinnedTreeStep() override = default;

    void describe() const override {
        RootDataAnalysisStep::describe();  // Вызываем базовую реализацию

        log("This step writes events into a 4D array of TTrees.", LogLevel::Info);
        log("Each TTree corresponds to a unique bin in (Q², W, cosθ, φ).", LogLevel::Info);
        log("Total number of trees: " + std::to_string(NUMBER_Q2 * NUMBER_W * NUMBER_COS_THETA * NUMBER_PHI), LogLevel::Debug);

        log("Stored variables per event:", LogLevel::Debug);
        log("  - MM (missing mass)", LogLevel::Debug);
        log("  - Q2 (momentum transfer squared)", LogLevel::Debug);
        log("  - W (invariant mass)", LogLevel::Debug);
        log("  - cos(theta) (cosine of pion angle in CM)", LogLevel::Debug);
        log("  - phi (azimuthal angle in CM)", LogLevel::Debug);

        log("Binning information:", LogLevel::Debug);
        log("  - Q2: " + std::to_string(NUMBER_Q2) + " bins (custom intervals)", LogLevel::Debug);
        log("  - W: " + std::to_string(NUMBER_W) + " bins from W_MIN=" + std::to_string(W_MIN) + " with step=" + std::to_string(STEP_W), LogLevel::Debug);
        log("  - cos(theta): " + std::to_string(NUMBER_COS_THETA) + " bins from -1 to 1, step=" + std::to_string(STEP_COS_THETA), LogLevel::Debug);
        log("  - phi: " + std::to_string(NUMBER_PHI) + " bins from 0 to 2π, step=" + std::to_string(STEP_PHI), LogLevel::Debug);
    }


protected:
    // --- Хранилище TTree'ов по всем ячейкам сетки ---
    std::vector<std::vector<std::vector<std::vector<TTree*>>>> trees;

    // --- Файл, в который будут сохранены выходные деревья ---
    // TFile* outputRootFile = nullptr;

    // --- Переменные, сохраняемые в деревья ---
    double MM = 0.0, Q2 = 0.0, W = 0.0, cos_theta = 0.0, phi_pi = 0.0;

    size_t skippedQ2W = 0;
    size_t skippedAngles = 0;

    size_t eventsProcessed = 0; ///< Сколько событий уже обработано
    size_t printEvery = 100000; ///< Частота логирования прогресса (по умолчанию)

    /**
     * @brief Инициализация: открытие файла, создание деревьев и связывание веток
     */
    bool initialize() override {
        if (!RootDataAnalysisStep::initialize()) {
            return false;
        }

        log("Initializing...", LogLevel::Info);

        // Создание структуры деревьев: 4D сетка по (Q², W, cosθ, φ)
        trees.resize(NUMBER_Q2);
        for (int i = 0; i < NUMBER_Q2; ++i) {
            trees[i].resize(NUMBER_W);
            for (int j = 0; j < NUMBER_W; ++j) {
                trees[i][j].resize(NUMBER_COS_THETA);
                for (int k = 0; k < NUMBER_COS_THETA; ++k) {
                    trees[i][j][k].resize(NUMBER_PHI);
                    for (int l = 0; l < NUMBER_PHI; ++l) {
                        std::string name = makeCellName(i, j, k, l);
                        std::string title = "TTree for bin Q2=" + std::to_string(i)
                                          + ", W=" + std::to_string(j)
                                          + ", CosTheta=" + std::to_string(k)
                                          + ", Phi=" + std::to_string(l);
                        trees[i][j][k][l] = new TTree(name.c_str(), title.c_str());
                        trees[i][j][k][l]->Branch("MM", &MM, "MM/D");
                        trees[i][j][k][l]->Branch("Q2", &Q2, "Q2/D");
                        trees[i][j][k][l]->Branch("W", &W, "W/D");
                        trees[i][j][k][l]->Branch("cos_theta", &cos_theta, "cos_theta/D");
                        trees[i][j][k][l]->Branch("phi", &phi_pi, "phi/D");
                        trees[i][j][k][l]->Branch("weight", &weight, "weight/D");
                    }
                }
            }
        }
        log("All output TTrees initialized (" + std::to_string(NUMBER_Q2 * NUMBER_W * NUMBER_COS_THETA * NUMBER_PHI) + " total).", LogLevel::Debug);

        return true;
    }

    /**
     * @brief Сохраняет все деревья в файл и закрывает его
     */
    bool finalize() override {
        log("Writing all TTrees to output file...", LogLevel::Info);
        outputFile->cd();
        for (int i = 0; i < NUMBER_Q2; ++i)
            for (int j = 0; j < NUMBER_W; ++j)
                for (int k = 0; k < NUMBER_COS_THETA; ++k)
                    for (int l = 0; l < NUMBER_PHI; ++l)
                        trees[i][j][k][l]->Write();

        log("Event processing summary:", LogLevel::Info);
        log("Total processed events: " + std::to_string(eventsProcessed) 
            + " out of " + std::to_string(totalEvents), LogLevel::Info);

        if (eventsProcessed > 0) {
            log("  Events skipped due to Q² or W out of range: " + std::to_string(skippedQ2W)
                + " (" + std::to_string(100.0 * skippedQ2W / eventsProcessed) + "%)", LogLevel::Info);

            log("  Events skipped due to angle bins (cosθ or φ): " + std::to_string(skippedAngles)
                + " (" + std::to_string(100.0 * skippedAngles / eventsProcessed) + "%)", LogLevel::Info);
        } else {
            log("No events were processed. Skipped stats are omitted.", LogLevel::Warning);
        }

        return RootDataAnalysisStep::finalize();
    }

    /**
     * @brief Обработка одного события: расчёт величин, биннинг, запись в нужное дерево
     */
    void processEvent() override {
        // Формируем 4-векторы электрона и пиона
        TLorentzVector electron, piPlus;
        electron.SetXYZM(e_px, e_py, e_pz, ELECTRON_MASS);
        piPlus.SetXYZM(pi_px, pi_py, pi_pz, CHARGED_PI_MASS);

        // Недостающий нейтрон и виртуальный фотон
        TLorentzVector neutron = BEAM + PROTON_TARGET - electron - piPlus;
        TLorentzVector gamma = BEAM - electron;

        // Расчёт переменных
        W = (BEAM + PROTON_TARGET - electron).M();
        Q2 = -gamma.M2();
        MM = neutron.M();

        // Определяем бин по Q² и W; пропускаем событие, если оно не попадает ни в один интервал
        int q2_bin = getQ2Bin(Q2);
        int w_bin = getWBin(W);
        if (q2_bin < 0 || w_bin < 0) {
            ++skippedQ2W;
            return;
        }

        // Переход в КС-систему
        TLorentzVector piPlusCM = piPlus;
        TVector3 betaLab(0.0, 0.0, -gamma.P() / (gamma.E() + PROTON_MASS));
        piPlusCM.RotateZ(-electron.Phi());
        piPlusCM.RotateY(gamma.Theta());
        piPlusCM.Boost(betaLab);
        piPlusCM.SetXYZM(-piPlusCM.Px(), -piPlusCM.Py(), -piPlusCM.Pz(), CHARGED_PI_MASS);

        cos_theta = TMath::Cos(piPlusCM.Theta());
        phi_pi = piPlusCM.Phi();

        // Биннинг по углам
        int cos_theta_bin = getCosThetaBin(cos_theta);
        int phi_bin = getPhiBin(phi_pi);
        if (cos_theta_bin < 0 || phi_bin < 0) {
            ++skippedAngles;
            return;
        }

        // Запись события в нужное дерево
        trees[q2_bin][w_bin][cos_theta_bin][phi_bin]->Fill();
    }

    void logProgress() override {
        ++eventsProcessed;
        if (eventsProcessed % printEvery == 0 || eventsProcessed == totalEvents) {
            double percent = 100.0 * eventsProcessed / totalEvents;
            log("Processed " + std::to_string(eventsProcessed) +
                " events (" + std::to_string(static_cast<int>(percent)) + "%)", LogLevel::Info);
        }
    }

private:
    /**
     * @brief Возвращает номер бина по W или -1, если вне диапазона
     */
    int getWBin(double w) {
        int bin = static_cast<int>((w - W_MIN) / STEP_W);
        return (bin < 0 || bin >= NUMBER_W) ? -1 : bin;
    }

    /**
     * @brief Возвращает номер бина по Q² из вручную заданных интервалов
     */
    int getQ2Bin(double q2) {
        for (int i = 0; i < NUMBER_Q2; ++i)
            if (q2 > STEPS_Q2[2 * i] && q2 < STEPS_Q2[2 * i + 1])
                return i;
        return -1;
    }

    /**
     * @brief Возвращает номер бина по cos(θ) или -1
     */
    int getCosThetaBin(double ct) {
        int bin = static_cast<int>((ct + 1.0) / STEP_COS_THETA);
        return (bin < 0 || bin >= NUMBER_COS_THETA) ? -1 : bin;
    }

    /**
     * @brief Возвращает номер бина по φ (в диапазоне 0–2π) или -1
     */
    int getPhiBin(double phi) {
        if (phi < 0.0) phi += 2 * M_PI;
        int bin = static_cast<int>(phi / STEP_PHI);
        return (bin < 0 || bin >= NUMBER_PHI) ? -1 : bin;
    }
};
