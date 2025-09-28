#include "hipoConversionStep.h"
#include "binnedTreeStep.h"
#include "histBuildefStepContrast.h"
#include "paramFitterStepExternal.h"

#include "drawHistStep.h"
#include "drawHistStepDebug.h"
#include "drawYieldStep.h"
#include "drawYieldSaplStep.h"

#include "mergeCyclesStep.h"
#include "source/ManageClasses/analysisManager.h"

//------------------------------------------------------------
// Основная конфигурация анализа
//------------------------------------------------------------
int main() {

    //--------------------------------------------------------
    // === Экспериментальные данные ===
    //--------------------------------------------------------
    auto stepConvertExp = std::make_shared<HipoConversionStep>(
        "convert_exp", "pass2.txt", "/volatile/clas12/borozdin/heavy_data/output.root"
    );

    auto stepBinExp = std::make_shared<BinnedTreeStep>(
        "bin_exp", "/volatile/clas12/borozdin/heavy_data/output.root", "/volatile/clas12/borozdin/heavy_data/binned.root"
    );

    auto stepHistExp = std::make_shared<HistBuilderStepContrast>(
        "hist_exp", "/volatile/clas12/borozdin/heavy_data/binned.root", "data_0808/histed.root"
    );

    auto stepFitExp = std::make_shared<ParamFitterStepExt>(
        "fit_exp", "data_0808/histed.root", "data_0808/fitted.root"
    );

    auto stepDrawHists = std::make_shared<DrawHistStep>(
        "draw_hists",
        std::map<std::string, std::string>{
            {"main", "data_0808/fitted.root"},
            {"hist", "data_0808/histed.root"}
        },
        "img/graphs_0808/"
    );

    auto stepDrawHistsDebug = std::make_shared<DrawHistStepDebug>(
        "draw_hists_debug", "data_0808/fitted.root", "data_0808/cond_debug.root"
    );

    auto stepDrawYield = std::make_shared<DrawYieldStep>(
        "draw_yield", "data_0808/fitted.root", "img/graphs_yield_0808/"
    );

    //--------------------------------------------------------
    // === Симуляционные данные ===
    //--------------------------------------------------------
    auto stepConvertSim = std::make_shared<HipoConversionStep>(
        "convert_sim", "pass3.txt", "/volatile/clas12/borozdin/heavy_data/outputSim.root", true
    );

    auto stepMergeSim = std::make_shared<MergeCyclesStep>(
        "merge_sim", "big_data/outputSim.root", "big_data/outputSimMerged.root"
    );

    auto stepBinSim = std::make_shared<BinnedTreeStep>(
        "bin_sim", "/volatile/clas12/borozdin/heavy_data/outputSim.root", "/volatile/clas12/borozdin/heavy_data/binnedSim.root"
    );

    auto stepHistSim = std::make_shared<HistBuilderStepContrast>(
        "hist_sim", "/volatile/clas12/borozdin/heavy_data/binnedSim.root", "data_0709/histedSim.root"
    );

    auto stepFitSim = std::make_shared<ParamFitterStepExt>(
        "fit_sim", "data_0709/histedSim.root", "data_0709/fittedSim.root"
    );

    auto stepDrawYieldExpVsSim = std::make_shared<DrawYieldSaplStep>(
        "draw_yield_exp_vs_sim",
        "data_0808/fitted.root",
        "data_0709/fittedSim.root",
        "img/graphs_yield_0709/"
    );

    //--------------------------------------------------------
    // === LUND данные ===
    //--------------------------------------------------------
    auto stepConvertLund = std::make_shared<HipoConversionLundStep>(
        "convert_lund", "pass3.txt", "/volatile/clas12/borozdin/heavy_data/outputLund.root", true
    );

    auto stepBinLund = std::make_shared<BinnedTreeStep>(
        "bin_lund", "/volatile/clas12/borozdin/heavy_data/outputLund.root", "/volatile/clas12/borozdin/heavy_data/binnedLund.root"
    );

    auto stepDrawYieldExpVsLund = std::make_shared<DrawYieldSaplExtendedStep>(
        "draw_yield_exp_vs_lund",
        "data_0808/fitted.root",
        "data_0709/fittedSim.root",
        "/volatile/clas12/borozdin/heavy_data/binnedLund.root",
        "img/graphs_yield_1809/"
    );

    //--------------------------------------------------------
    // === Рекалькуляции и домашние данные ===
    //--------------------------------------------------------
    auto stepBinHome = std::make_shared<BinnedTreeStep>(
        "bin_home", "big_data/output.root", "big_data/binned_home.root"
    );

    auto stepHistHome = std::make_shared<HistBuilderStepContrast>(
        "hist_home", "big_data/binned.root", "big_data/histed_home.root"
    );

    auto stepFitExpRecalc = std::make_shared<ParamFitterStepExt>(
        "fit_exp_recalc", "data_0808/histed.root", "data_0808/fitted.root"
    );

    auto stepFitSimRecalc = std::make_shared<ParamFitterStepExt>(
        "fit_sim_recalc", "data_0709/histedSim.root", "data_0709/fittedSim.root"
    );

    //--------------------------------------------------------
    // === Ветви анализа ===
    //--------------------------------------------------------
    auto branchExp = std::make_shared<AnalysisBranch>();
    // branchExp->addStep(stepConvertExp);
    branchExp->addStep(stepBinExp);
    branchExp->addStep(stepHistExp);
    branchExp->addStep(stepFitExp);
    // branchExp->addStep(stepDrawHists);
    // branchExp->addStep(stepDrawYield);

    auto branchSim = std::make_shared<AnalysisBranch>();
    branchSim->addStep(stepConvertSim);
    branchSim->addStep(stepBinSim);
    branchSim->addStep(stepHistSim);
    branchSim->addStep(stepFitSim);
    // branchSimYield->addStep(stepDrawYieldExpVsSim);

    auto branchSimMerged = std::make_shared<AnalysisBranch>();
    branchSimMerged->addStep(stepMergeSim);

    auto branchLund = std::make_shared<AnalysisBranch>();
    branchLund->addStep(stepConvertLund);
    branchLund->addStep(stepBinLund);
    branchLund->addStep(stepDrawYieldExpVsLund);

    auto branchRecalc = std::make_shared<AnalysisBranch>();
    branchRecalc->addStep(stepFitExpRecalc);
    branchRecalc->addStep(stepFitSimRecalc);

    //--------------------------------------------------------
    // === Менеджер анализа ===
    //--------------------------------------------------------
    AnalysisManager manager;
    manager.addBranch("experiment", branchExp);
    manager.addBranch("simulation_raw", branchSim);
    manager.addBranch("simulation_merge", branchSimMerged);
    manager.addBranch("lund", branchLund);
    manager.addBranch("recalc", branchRecalc);

    manager.describe();

    manager.runBranch("experiment");
    manager.runBranch("simulation_raw");
    // manager.runBranch("simulation_merge");
    manager.runBranch("lund");
    // manager.runBranch("recalc");

    return 0;
}
