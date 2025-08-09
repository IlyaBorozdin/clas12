#include "hipoConversionStep.h"
#include "binnedTreeStep.h"
#include "histBuildefStepContrast.h"
#include "paramFitterStepExternal.h"
#include "drawHistStep.h"
#include "drawHistStepDebug.h"
#include "drawYieldStep.h"
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

    AnalysisManager manager;
    manager.addBranch("experimental", branch2);
    manager.addBranch("yield builder", branch3);
    manager.addBranch("output simulation", branch4);
    manager.addBranch("home", branch5);

    manager.describe();
    manager.runBranch("experimental");
    // manager.runBranch("yield builder");
    // manager.runBranch("output simulation");
    // manager.runBranch("home");

    return 0;
}