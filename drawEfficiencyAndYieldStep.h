#pragma once

#include "drawYieldSaplStep.h"
#include "utils/fittingYield.h"

class DrawEfficiencyAndYieldStep : public DrawYieldSaplExtendedStep {
public:
    DrawEfficiencyAndYieldStep(const std::string& stepName,
                               const std::string& inputFileName,
                               const std::string& inputFileSimName,
                               const std::string& inputFileLundName,
                               const std::string& outputDirName = "img/")
        : DrawYieldSaplExtendedStep(stepName, inputFileName, inputFileSimName, inputFileLundName, outputDirName) {}
    
protected:
    std::vector<std::vector<TCanvas*>> effCanvases;
    std::vector<std::vector<std::vector<TGraphErrors*>>> graphsEff;
    std::vector<std::vector<std::vector<TGraphErrors*>>> graphsCorr;

    bool initialize() override {
        if (!DrawYieldSaplStep::initialize()) {
            return false;
        }

        effCanvases.resize(NUMBER_Q2);
        graphsEff.resize(NUMBER_Q2);
        graphsCorr.resize(NUMBER_Q2);
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            effCanvases[iq2].resize(NUMBER_W, nullptr);
            graphsEff[iq2].resize(NUMBER_W);
            graphsCorr[iq2].resize(NUMBER_W);
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                std::string cname = "eff_canvas_q2_" + std::to_string(iq2 + 1) + "_w_" + std::to_string(iw + 1);
                effCanvases[iq2][iw] = new TCanvas(cname.c_str(), cname.c_str(), 2000, 800);
                effCanvases[iq2][iw]->Divide(5, 2);

                graphsEff[iq2][iw].resize(NUMBER_COS_THETA, nullptr);
                graphsCorr[iq2][iw].resize(NUMBER_COS_THETA, nullptr);
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    graphsEff[iq2][iw][ict] = new TGraphErrors();
                    graphsEff[iq2][iw][ict]->SetTitle(generateTitleForEfficiencyPlot(iq2, iw, ict).c_str());
                    graphsEff[iq2][iw][ict]->SetLineColor(kBlue);
                    graphsEff[iq2][iw][ict]->SetLineWidth(2);

                    graphsCorr[iq2][iw][ict] = new TGraphErrors();
                    graphsCorr[iq2][iw][ict]->SetTitle(generateTitleForYieldPlot(iq2, iw, ict).c_str());
                    graphsCorr[iq2][iw][ict]->SetLineColor(kBlue);
                    graphsCorr[iq2][iw][ict]->SetLineWidth(2);
                }
            }
        }

        return true;
    }

    bool run() override {
        experiment.safeRun();
        simulation.safeRun();
        lund.safeRun();

        graphsLund = lund.getGraphs();
        graphsExp = experiment.getGraphs();
        graphsSim = simulation.getGraphs();

        return (lund.getStatus() == StepStatus::Done || lund.getStatus() == StepStatus::Skipped) &&
               (experiment.getStatus() == StepStatus::Done || experiment.getStatus() == StepStatus::Skipped) &&
               (simulation.getStatus() == StepStatus::Done || simulation.getStatus() == StepStatus::Skipped);
    }

    bool finalize() override {
        for (int iq2 = 0; iq2 < NUMBER_Q2; ++iq2) {
            for (int iw = 0; iw < NUMBER_W; ++iw) {
                for (int ict = 0; ict < NUMBER_COS_THETA; ++ict) {
                    TGraphErrors* grExp = graphsExp[iq2][iw][ict];  // эксперимент
                    TGraphErrors* grCorr = graphsCorr[iq2][iw][ict];  // эксперимент
                    TGraphErrors* grSim = graphsSim[iq2][iw][ict];  // симуляция
                    TGraphErrors* grLund = graphsLund[iq2][iw][ict];// LUND
                    TGraphErrors* grEff = graphsEff[iq2][iw][ict];  // эффективность

                    for (int ip = 0; ip < NUMBER_PHI; ++ip) {
                        double phi = (ip + 0.5) * (2 * M_PI / NUMBER_PHI);

                        int idxExp  = findPointIndexByPhi(grExp, phi);
                        int idxSim = findPointIndexByPhi(grSim, phi);
                        int idxLund = findPointIndexByPhi(grLund, phi);
                        if (idxLund < 0 || idxSim < 0) continue;

                        double phiLund, yLund;
                        grLund->GetPoint(idxLund, phiLund, yLund);
                        double errLund = grLund->GetErrorY(idxLund);

                        double phiSim, ySim;
                        grSim->GetPoint(idxSim, phiSim, ySim);
                        double errSim = grSim->GetErrorY(idxSim);

                        if (yLund > 0.0 && ySim > 0.0) {
                            double eff = ySim / yLund;
                            double errEff = eff * std::sqrt(std::pow(errLund / yLund, 2) + std::pow(errSim / ySim, 2));

                            int nEff = grEff->GetN();
                            grEff->SetPoint(nEff, phi, eff);
                            grEff->SetPointError(nEff, M_PI / NUMBER_PHI, errEff);

                            if (idxExp < 0) continue;
                            
                            double phiExp, yExp;
                            grExp->GetPoint(idxExp, phiExp, yExp);
                            double errExp = grExp->GetErrorY(idxExp);

                            if (yExp > 0) {
                                double yCorr = yExp / eff;
                                double errCorr = yCorr * std::sqrt(std::pow(errEff / eff, 2) + std::pow(errExp / yExp, 2));

                                int nCorr = grCorr->GetN();
                                grCorr->SetPoint(nCorr, phi, yCorr);
                                grCorr->SetPointError(nCorr, M_PI / NUMBER_PHI, errCorr);
                            }
                        }
                    }

                    // отрисовка на канвасе
                    effCanvases[iq2][iw]->cd(ict+1);
                    grEff->Draw("APL");

                    TF1* fitFunc = fitting.fit(grCorr, iq2, iw, ict);
                    fitFunc->SetLineColor(kBlue);

                    canvases[iq2][iw]->cd(ict + 1);
                    grCorr->Draw("AP");
                    fitFunc->Draw("same");
                    gPad->Update(); // важно: создаёт скрытую гистограмму осей

                    auto* hist = grCorr->GetHistogram();
                    if (hist) {
                        hist->SetMinimum(0); // нижняя граница оси Y = 0
                    }

                    delete fitFunc;

                }

                effCanvases[iq2][iw]->Update();
                std::string effFigName = "eff_q2_" + std::to_string(iq2 + 1) +
                                         "_w_" + std::to_string(iw + 1) + ".png";
                effCanvases[iq2][iw]->SaveAs((outputDirName + effFigName).c_str());

                canvases[iq2][iw]->Update();
                std::string figName = "yield_q2_" + std::to_string(iq2 + 1) +
                                      "_w_" + std::to_string(iw + 1) + ".png";
                canvases[iq2][iw]->SaveAs((outputDirName + figName).c_str());
            }
        }

        return true;
    }

private:
    FittingYield fitting;

    int findPointIndexByPhi(TGraphErrors* gr, double phi, double tolerance = 1e-6) {
        int nPoints = gr->GetN();
        for (int i = 0; i < nPoints; ++i) {
            double x, y;
            gr->GetPoint(i, x, y);
            if (std::fabs(x - phi) < tolerance) {
                return i; // нашли точку
            }
        }
        return -1; // не нашли
    }

};