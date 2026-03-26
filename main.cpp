#include "hipoConversionStep.h"
#include "binnedTreeStep.h"
#include "histBuilderStepContrast.h"
#include "paramFitterStepExternal.h"

#include "drawHistStep.h"
#include "drawHistStepDebug.h"
#include "drawYieldStep.h"
#include "drawYieldSaplStep.h"
#include "drawEfficiencyAndYieldStep.h"

#include "mergeCyclesStep.h"
#include "mergeFilesStep.h"
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
        "bin_exp", std::string(HEAVY_BASE_PATH) + "/output.root", std::string(HEAVY_BASE_PATH) + "/binned_1703.root"
    );

    auto stepHistExp = std::make_shared<HistBuilderStepContrast>(
        "hist_exp", std::string(HEAVY_BASE_PATH) + "/binned_1703.root", "data_1703/histed.root"
    );

    auto stepFitExp = std::make_shared<ParamFitterStepExt>(
        "fit_exp", "data_1703/histed.root", "data_1703/fitted.root", "/home/borozdin/move_to_ifarm_1903_1/utils/fit_configs.json"
    );

    auto stepDrawHists = std::make_shared<DrawHistStep>(
        "draw_hists",
        std::map<std::string, std::string>{
            {"main", "data_1703/fitted.root"},
            {"hist", "data_1703/histed.root"}
        },
        "img/graphs_1703/exp/"
    );

    auto stepDrawYield = std::make_shared<DrawYieldStep>(
        "draw_yield", "data_0808/fitted.root", "img/graphs_yield_0808/"
    );

    //--------------------------------------------------------
    // === Симуляционные данные ===
    //--------------------------------------------------------
    auto stepConvertSim = std::make_shared<HipoConversionSimStep>(
        "convert_sim", "pass3.txt", std::string(HEAVY_BASE_PATH) + "/outputSim_0611.root", true
    );

    auto stepMergeSim = std::make_shared<MergeCyclesStep>(
        "merge_sim", "/volatile/clas12/borozdin/heavy_data/outputSim.root", "/volatile/clas12/borozdin/heavy_data/outputSimC.root"
    );

    auto stepMergeFileSim = std::make_shared<MergeFilesStep>(
        "merge_file_sim",
        std::map<std::string, std::string>{
            {"main", std::string(HEAVY_BASE_PATH) + "/outputSim_0611.root"},
            {"add1", std::string(HEAVY_BASE_PATH) + "/outputSim_0711.root"}
        },
        // "/volatile/clas12/borozdin/heavy_data/outputSimMerged_3110.root"
        std::string(HEAVY_BASE_PATH) + "/outputSimMerged_0711.root"
    );

    auto stepBinSim = std::make_shared<BinnedTreeStep>(
        "bin_sim", std::string(HEAVY_BASE_PATH) + "/outputSimMerged_0711.root", std::string(HEAVY_BASE_PATH) + "/binnedSim_1703.root"
    );

    auto stepHistSim = std::make_shared<HistBuilderStepContrast>(
        "hist_sim", std::string(HEAVY_BASE_PATH) + "/binnedSim_1703.root", "data_1703/histedSim.root"
    );

    auto stepFitSim = std::make_shared<ParamFitterStepExt>(
        "fit_sim", "data_1703/histedSim.root", "data_1703/fittedSim.root", "/home/borozdin/move_to_ifarm_1903_1/utils/fit_configs.json"
    );

    auto stepDrawHistsSim = std::make_shared<DrawHistStep>(
        "draw_hists_sim",
        std::map<std::string, std::string>{
            {"main", "data_1703/fittedSim.root"},
            {"hist", "data_1703/histedSim.root"}
        },
        "img/graphs_1703/sim/"
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
        "convert_lund", "pass3.txt", std::string(HEAVY_BASE_PATH) + "/outputLund_0711.root", true
    );

    auto stepMergeLund = std::make_shared<MergeCyclesStep>(
        "merge_lund", "/volatile/clas12/borozdin/heavy_data/outputLundB.root", "/volatile/clas12/borozdin/heavy_data/outputLundC.root"
    );

    auto stepMergeFileLund = std::make_shared<MergeFilesStep>(
        "merge_file_lund",
        std::map<std::string, std::string>{
            {"main", std::string(HEAVY_BASE_PATH) + "/outputLund_0611.root"},
            {"add1", std::string(HEAVY_BASE_PATH) + "/outputLund_0711.root"}
        },
        std::string(HEAVY_BASE_PATH) + "/outputLundMerged_0711.root"
    );

    auto stepBinLund = std::make_shared<BinnedTreeStep>(
        "bin_lund", std::string(HEAVY_BASE_PATH) + "/outputLundMerged_0711.root", std::string(HEAVY_BASE_PATH) + "/binnedLund_1703.root"
    );

    auto stepHistLund = std::make_shared<HistBuilderStepContrast>(
        "hist_lund", std::string(HEAVY_BASE_PATH) + "/binnedLund_1703.root", "data_1703/histedLund.root"
    );

    auto stepDrawYieldExpVsSimVsLund = std::make_shared<DrawYieldSaplExtendedStep>(
        "draw_yield_exp_vs_sim_vs_lund",
        "data_0808/fitted.root",
        "data_0709/fittedSim.root",
        "data_0709/histedLund.root",
        "img/graphs_yield_vs_1809/"
    );

    auto stepDrawEfficiencyAndYieldExpVsLund = std::make_shared<DrawEfficiencyAndYieldStep>(
        "draw_yield_exp_vs_lund",
        "data_1703/fitted.root",
        "data_1703/fittedSim.root",
        "data_1703/histedLund.root",
        "img/graphs_yield_1703/"
    );

    //--------------------------------------------------------
    // === Ветви анализа ===
    //--------------------------------------------------------
    auto branchExp = std::make_shared<AnalysisBranch>();
    // branchExp->addStep(stepConvertExp);
    branchExp->addStep(stepBinExp);
    branchExp->addStep(stepHistExp);
    branchExp->addStep(stepFitExp);
    branchExp->addStep(stepDrawHists);
    // branchExp->addStep(stepDrawYield);

    auto branchSim = std::make_shared<AnalysisBranch>();
    // branchSim->addStep(stepConvertSim);
    branchSim->addStep(stepBinSim);
    branchSim->addStep(stepHistSim);
    branchSim->addStep(stepFitSim);
    branchSim->addStep(stepDrawHistsSim);
    // branchSim->addStep(stepDrawYieldExpVsSim);

    auto branchLund = std::make_shared<AnalysisBranch>();
    // branchLund->addStep(stepConvertLund);
    branchLund->addStep(stepBinLund);
    branchLund->addStep(stepHistLund);
    // branchLund->addStep(stepDrawYieldExpVsSimVsLund);
    branchLund->addStep(stepDrawEfficiencyAndYieldExpVsLund);

    auto branchRecord = std::make_shared<AnalysisBranch>();
    // branchRecord->addStep(stepConvertSim);
    branchRecord->addStep(stepConvertLund);

    auto branchMerge = std::make_shared<AnalysisBranch>();
    // branchMerge->addStep(stepMergeSim);
    // branchMerge->addStep(stepMergeLund);
    branchMerge->addStep(stepMergeFileSim);
    branchMerge->addStep(stepMergeFileLund);

    //--------------------------------------------------------
    // === Менеджер анализа ===
    //--------------------------------------------------------
    AnalysisManager manager;
    manager.addBranch("experiment", branchExp);
    manager.addBranch("simulation", branchSim);
    manager.addBranch("lund", branchLund);
    manager.addBranch("record", branchRecord);
    manager.addBranch("merge", branchMerge);

    manager.describe();

    // manager.runBranch("record");
    // manager.runBranch("merge");

    manager.runBranch("experiment");
    manager.runBranch("simulation");
    manager.runBranch("lund");

    return 0;
}
