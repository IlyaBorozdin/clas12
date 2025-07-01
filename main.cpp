#include "hipoConversionStep.h"
#include "binnedTreeStep.h"
#include "histBuilderStep.h"
#include "paramFitterStep.h"
#include "drawHistStep.h"
#include "source/ManageClasses/analysisManager.h"

int main() {

    auto A = std::make_shared<HipoConversionStep>(
        "A", "iv/pass2.txt", "data/output.root"
    );
    auto B = std::make_shared<BinnedTreeStep>(
        "B", "data/output.root", "data/binned.root"
    );
    auto C = std::make_shared<HistBuilderStep>(
        "C", "data/binned.root", "data/histed.root"
    );
    auto D = std::make_shared<ParamFitterStep>(
        "D", "data/histed.root", "data/fit_data.root"
    );
    auto E = std::make_shared<DrawHistStep>(
        "E", std::map<std::string, std::string>{ {"main", "data/fit_data.root"}, {"hist", "data/histed.root"} }, "img/graphs/"
    );

    auto branch = std::make_shared<AnalysisBranch>();
    branch->addStep(A);
    branch->addStep(B);
    branch->addStep(C);
    branch->addStep(D);
    branch->addStep(E);

    AnalysisManager manager;
    manager.addBranch("default", branch);

    manager.describe();
    manager.runBranch("default");

    return 0;
}