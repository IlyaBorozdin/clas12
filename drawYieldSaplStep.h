#pragma once

#include "drawYieldStep.h"
#include "drawYieldLundShortStep.h"

class DrawYieldShortStep : public DrawYieldStep {
public:
    DrawYieldShortStep(const std::string& stepName,
                       const std::string& inputFileName)
        : DrawYieldStep(stepName, inputFileName) {}

    virtual ~DrawYieldShortStep() = default;

    std::vector<std::vector<std::vector<TGraphErrors*>>> getGraphs() const {
        return graphs;
    }

protected:
    bool finalize() override {
        log("Yield plotting complete.", LogLevel::Info);
        log("Processed entries: " + std::to_string(eventsProcessed), LogLevel::Info);

        return RootFitAnalysisStep::finalize();
    }
};

class DrawYieldSaplStep : public BaseStep {
public:
    DrawYieldSaplStep(const std::string& stepName,
                      const std::string& inputFileName,
                      const std::string& inputFileSimName,
                      const std::string& outputDirName = "img/")
        : BaseStep(stepName),
          experiment(stepName, inputFileName),
          simulation(stepName, inputFileSimName),
          outputDirName(outputDirName) {}

protected:
    std::string outputDirName;

    DrawYieldShortStep experiment;
    DrawYieldShortStep simulation;

    std::vector<std::vector<std::vector<TGraphErrors*>>> graphsExp;
    std::vector<std::vector<std::vector<TGraphErrors*>>> graphsSim;

    std::vector<std::vector<TCanvas*>> canvases;

    bool initialize() override {
        canvases.resize(NUMBER_Q2);
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            canvases[iq2].resize(NUMBER_W, nullptr);
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                std::string cname = "snd_canvas_q2_" + std::to_string(iq2 + 1) + "_w_" + std::to_string(iw + 1);
                canvases[iq2][iw] = new TCanvas(cname.c_str(), cname.c_str(), 1600, 1200);
                canvases[iq2][iw]->Divide(4, 3);
            }
        }

        return true;
    }

    bool run() override {
        experiment.safeRun();
        simulation.safeRun();

        graphsExp = experiment.getGraphs();
        graphsSim = simulation.getGraphs();

        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    graphsSim[iq2][iw][ict]->SetLineColor(kRed);
                }
            }
        }

        return (experiment.getStatus() == StepStatus::Done || experiment.getStatus() == StepStatus::Skipped) &&
               (simulation.getStatus() == StepStatus::Done || simulation.getStatus() == StepStatus::Skipped);
    }

    bool finalize() override {
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    canvases[iq2][iw]->cd(ict + 1);
                    graphsExp[iq2][iw][ict]->Draw("APL");
                    graphsSim[iq2][iw][ict]->Draw("PL SAME");
                }

                canvases[iq2][iw]->Update();

                std::string figName = "yield_q2_" + std::to_string(iq2 + 1) +
                                    "_w_" + std::to_string(iw + 1) + ".png";
                canvases[iq2][iw]->SaveAs((outputDirName + figName).c_str());
            }
        }
        return true;
    }
};

class DrawYieldSaplExtendedStep : public DrawYieldSaplStep {
public:
    DrawYieldSaplExtendedStep(const std::string& stepName,
                              const std::string& inputFileName,
                              const std::string& inputFileSimName,
                              const std::string& inputFileLundName,
                              const std::string& outputDirName = "img/")
        : DrawYieldSaplStep(stepName, inputFileName, inputFileSimName, outputDirName),
          lund(stepName, inputFileLundName) {}

protected:
    DrawYieldLundShortStep lund;
    std::vector<std::vector<std::vector<TGraphErrors*>>> graphsLund;

    bool run() override {
        // запускаем базовый run() (эксперимент + симуляция)
        if (!DrawYieldSaplStep::run()) return false;

        // запускаем LUND
        lund.safeRun();
        graphsLund = lund.getGraphs();

        // нормировка графиков
        normalizeGraphs(graphsExp);
        normalizeGraphs(graphsSim);
        normalizeGraphs(graphsLund);

        return (lund.getStatus() == StepStatus::Done || lund.getStatus() == StepStatus::Skipped);
    }

    bool finalize() override {
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    canvases[iq2][iw]->cd(ict + 1);
                    graphsExp[iq2][iw][ict]->Draw("APL");
                    graphsSim[iq2][iw][ict]->Draw("PL SAME");
                    graphsLund[iq2][iw][ict]->Draw("PL SAME");
                }

                canvases[iq2][iw]->Update();

                std::string figName = "yield_q2_" + std::to_string(iq2 + 1) +
                                    "_w_" + std::to_string(iw + 1) + ".png";
                canvases[iq2][iw]->SaveAs((outputDirName + figName).c_str());
            }
        }
        return true;
    }

private:
    void normalizeGraphs(std::vector<std::vector<std::vector<TGraphErrors*>>>& graphs) {
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    TGraphErrors* gr = graphs[iq2][iw][ict];
                    if (!gr) continue;

                    double sum = 0.0;
                    int n = gr->GetN();
                    for (int i = 0; i < n; ++i) {
                        double x, y;
                        gr->GetPoint(i, x, y);
                        sum += y;
                    }

                    if (sum > 0) {
                        for (int i = 0; i < n; ++i) {
                            double x, y;
                            gr->GetPoint(i, x, y);
                            gr->SetPoint(i, x, y / sum);
                            gr->SetPointError(i, gr->GetErrorX(i), gr->GetErrorY(i) / sum);
                        }
                    }
                }
            }
        }
    }
};
