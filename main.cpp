#include "hipoConversionStep.h"
#include "binnedTreeStep.h"
#include "histBuilderStep.h"
#include "histBuildefStepContrast.h"
#include "paramFitterStep.h"
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
        "B", "output.root", "data_0307/binned.root"
    );

    auto C = std::make_shared<HistBuilderStep>(
        "C", "data/binned.root", "data/histed.root"
    );
    auto F = std::make_shared<HistBuilderStepContrast>(
        // "F", "data/binned.root", "data/histed.root"
        "F", "data_0307/binned.root", "data_0307/histed.root"
    );

    auto D = std::make_shared<ParamFitterStep>(
        "D", "data/histed.root", "data/fit_data.root"
    );
    auto H = std::make_shared<ParamFitterStepExt>(
        // "H", "data/histed.root", "data/fit_data.root"
        "H", "data_0307/histed.root", "data_0307/fitted.root"
    );

    auto E = std::make_shared<DrawHistStep>(
        // "E", std::map<std::string, std::string>{ {"main", "data/fit_data.root"}, {"hist", "data/histed.root"} }, "img/graphs/"
        "E", std::map<std::string, std::string>{ {"main", "data_0307/fitted.root"}, {"hist", "data_0307/histed.root"} }, "img/graphs_1807/"
    );
    auto I = std::make_shared<DrawHistStepDebug>(
        // "I", "data/fit_data.root", "data/cond_debug.root"
        "I", "data_0307/fitted.root", "data_0307/cond_debug.root"
    );

    auto G = std::make_shared<DrawYieldStep>(
        // "E", std::map<std::string, std::string>{ {"main", "data/fit_data.root"}, {"hist", "data/histed.root"} }, "img/graphs/"
        "G", "data_0307/fitted.root", "img/graphs_yield_1807/"
    );
    

    auto branch1 = std::make_shared<AnalysisBranch>();
    branch1->addStep(A);
    branch1->addStep(B);
    branch1->addStep(C);
    branch1->addStep(D);
    branch1->addStep(E);

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

    AnalysisManager manager;
    manager.addBranch("default", branch1);
    manager.addBranch("experimental", branch2);
    manager.addBranch("yield builder", branch3);
    manager.addBranch("output simulation", branch4);

    manager.describe();
    // manager.runBranch("default");
    // manager.runBranch("experimental");
    // manager.runBranch("yield builder");
    manager.runBranch("output simulation");

    return 0;
}