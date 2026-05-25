#pragma once

#include "TStyle.h"

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

        gStyle->SetTitleFont(42, "t");
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTextFont(42);
        gStyle->SetTitleSize(0.05, "XYZ"); // Сделаем чуть крупнее для TGraph
        gStyle->SetLabelSize(0.04, "XYZ");

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
                    effCanvases[iq2][iw]->GetPad(ict + 1)->SetGrid();

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
                    // 1. Настройка и отрисовка эффективности
                    // --- 1. Отрисовка эффективности ---
                    effCanvases[iq2][iw]->cd(ict + 1);
                    grEff->Draw("APL"); 
                    gPad->Update(); // Создаем оси

                    double minX_e, maxX_e, minY_e, maxY_e;
                    grEff->ComputeRange(minX_e, minY_e, maxX_e, maxY_e);
                    if (grEff->GetHistogram()) {
                        grEff->GetHistogram()->SetMaximum(maxY_e * 1.35);
                        grEff->GetHistogram()->SetMinimum(0);
                    }
                    addKinematicInfoToGraph(grEff, iq2, iw, ict);
                    gPad->Modified(); // Помечаем, что пад изменился

                    // --- 2. Отрисовка выхода (Yield) ---
                    canvases[iq2][iw]->cd(ict + 1);
                    grCorr->Draw("AP");

                    TF1* fitFunc = fitting.fit(grCorr, iq2, iw, ict);
                    fitFunc->SetLineColor(kRed);
                    fitFunc->Draw("same");

                    gPad->Update(); 

                    double minX_c, maxX_c, minY_c, maxY_c;
                    grCorr->ComputeRange(minX_c, minY_c, maxX_c, maxY_c);
                    if (grCorr->GetHistogram()) {
                        grCorr->GetHistogram()->SetMaximum(maxY_c * 1.35);
                        grCorr->GetHistogram()->SetMinimum(0);
                    }
                    addKinematicInfoToGraph(grCorr, iq2, iw, ict);
                    gPad->Modified();

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

    void addKinematicInfoToGraph(TGraphErrors* gr, int iq2, int iw, int ict) const {
        // Рассчитываем границы (аналогично вашему HistBuilderStep)
        double q2Min = STEPS_Q2[2 * iq2], q2Max = STEPS_Q2[2 * iq2 + 1];
        double wMin = W_MIN + STEP_W * iw, wMax = W_MIN + STEP_W * (iw + 1);
        double cosMin = -1.0 + STEP_COS_THETA * ict, cosMax = -1.0 + STEP_COS_THETA * (ict + 1);

        // Позиционируем в верхнем левом углу, чтобы не перекрывать точки (обычно они справа)
        TPaveText* pt = new TPaveText(0.15, 0.70, 0.45, 0.88, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(0);
        pt->SetTextAlign(12);
        pt->SetTextFont(42);
        pt->SetTextSize(0.035);

        pt->AddText(Form("Q^{2}: %.2f-%.2f GeV^{2}", q2Min, q2Max));
        pt->AddText(Form("W: %.3f-%.3f GeV", wMin, wMax));
        pt->AddText(Form("cos(#theta^{*}): %.2f-%.2f", cosMin, cosMax));

        // Добавляем объект в список функций графа, чтобы он рисовался вместе с ним
        gr->GetListOfFunctions()->Add(pt);
    }

};