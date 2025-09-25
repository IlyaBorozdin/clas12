#include "hipoConversionStep.h"
#include "binnedTreeStep.h"
#include "histBuildefStepContrast.h"
#include "paramFitterStepExternal.h"
#include "drawHistStep.h"
#include "drawHistStepDebug.h"
#include "drawYieldStep.h"
#include "mergeCyclesStep.h"
#include "drawYieldSaplStep.h"
#include "source/ManageClasses/analysisManager.h"

int main() {

    auto A = std::make_shared<HipoConversionStep>(
        // "A", "iv/pass2.txt", "data/output.root"
        "A", "pass2.txt", "output.root"
    );
    auto J = std::make_shared<HipoConversionStep>(
        "J", "pass3.txt", "outputSim.root", true
    );

    auto B = std::make_shared<BinnedTreeStep>(
        // "B", "data/output.root", "data/binned.root"
        "B", "output.root", "data_0808/binned.root"
    );
    auto B_home = std::make_shared<BinnedTreeStep>(
        // "B", "data/output.root", "data/binned.root"
        "B", "big_data/output.root", "big_data/binned_home.root"
    );

    auto F = std::make_shared<HistBuilderStepContrast>(
        // "F", "data/binned.root", "data/histed.root"
        "F", "data_0808/binned.root", "data_0808/histed.root"
    );
    auto F_home = std::make_shared<HistBuilderStepContrast>(
        // "F", "data/binned.root", "data/histed.root"
        "F", "big_data/binned.root", "big_data/histed_home.root"
    );

    auto H = std::make_shared<ParamFitterStepExt>(
        // "H", "data/histed.root", "data/fit_data.root"
        "H", "data_0808/histed.root", "data_0808/fitted.root"
    );

    auto E = std::make_shared<DrawHistStep>(
        // "E", std::map<std::string, std::string>{ {"main", "data/fit_data.root"}, {"hist", "data/histed.root"} }, "img/graphs/"
        "E", std::map<std::string, std::string>{ {"main", "data_0808/fitted.root"}, {"hist", "data_0808/histed.root"} }, "img/graphs_0808/"
    );
    auto I = std::make_shared<DrawHistStepDebug>(
        // "I", "data/fit_data.root", "data/cond_debug.root"
        "I", "data_0808/fitted.root", "data_0808/cond_debug.root"
    );

    auto G = std::make_shared<DrawYieldStep>(
        // "E", std::map<std::string, std::string>{ {"main", "data/fit_data.root"}, {"hist", "data/histed.root"} }, "img/graphs/"
        "G", "data_0808/fitted.root", "img/graphs_yield_0808/"
    );

    auto A_SIM = std::make_shared<MergeCyclesStep>(
        "A_SIM", "big_data/outputSim.root", "big_data/outputSimMerged.root"
    );

    auto B_SIM = std::make_shared<BinnedTreeStep>(
        "B_SIM", "big_data/outputSimMerged.root", "big_data/binnedSim.root"
    );

    auto F_SIM = std::make_shared<HistBuilderStepContrast>(
        "F_SIM", "data_0709/binnedSim.root", "data_0709/histedSim.root"
    );

    auto H_SIM = std::make_shared<ParamFitterStepExt>(
        "H", "data_0709/histedSim.root", "data_0709/fittedSim.root"
    );

    auto G_SIM = std::make_shared<DrawYieldSaplStep>(
        "G_SIN", "data_0808/fitted.root", "data_0709/fittedSim.root", "img/graphs_yield_0709/"
    );

    auto J_LUND = std::make_shared<HipoConversionLundStep>(
        "J_LUND", "pass4.txt", "outputLund.root", true
    );

    auto B_LUND = std::make_shared<BinnedTreeStep>(
        "B_LUND", "outputLund.root", "data_1809/binned.root"
    );

    auto G_LUND = std::make_shared<DrawYieldSaplExtendedStep>(
        "G_LUND", "data_0808/fitted.root", "data_0709/fittedSim.root", "data_1809/binned.root", "img/graphs_yield_1809/"
    );


    auto branch2 = std::make_shared<AnalysisBranch>();
    branch2->addStep(A);
    branch2->addStep(B);
    // С этого момента изменения !!! (F и H)
    branch2->addStep(F);
    branch2->addStep(H);
    branch2->addStep(E);
    // branch2->addStep(I);
    branch2->addStep(G);

    auto branch3 = std::make_shared<AnalysisBranch>();
    branch3->addStep(G);

    auto branch4 = std::make_shared<AnalysisBranch>();
    branch4->addStep(J);

    auto branch5 = std::make_shared<AnalysisBranch>();
    branch5->addStep(F_home);

    auto branchSIM6 = std::make_shared<AnalysisBranch>();
    // branchSIM6->addStep(B_SIM);
    // branchSIM6->addStep(F_SIM);
    // branchSIM6->addStep(H_SIM);
    branchSIM6->addStep(G_SIM);

    auto branchSIM7 = std::make_shared<AnalysisBranch>();
    branchSIM7->addStep(A_SIM);

    auto branchLUND8 = std::make_shared<AnalysisBranch>();
    branchLUND8->addStep(J_LUND);
    // branchLUND8->addStep(B_LUND);
    // branchLUND8->addStep(G_LUND);

    AnalysisManager manager;
    manager.addBranch("experimental", branch2);
    manager.addBranch("yield builder", branch3);
    manager.addBranch("output simulation", branch4);
    manager.addBranch("home", branch5);

    manager.addBranch("home sim", branchSIM6);
    manager.addBranch("merge sim", branchSIM7);
    manager.addBranch("yield lund", branchLUND8);

    manager.describe();
    // manager.runBranch("experimental");
    // manager.runBranch("yield builder");
    // manager.runBranch("output simulation");
    // manager.runBranch("home");
    // manager.runBranch("home sim");
    manager.runBranch("yield lund");

    return 0;
}