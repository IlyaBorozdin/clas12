#include <string>
#include <iostream>
#include <fstream>
#include <vector>

//#include "neutron_piPlus.h"
//#include "yield_in_q2-w.h"
//#include "choiceRoot.h"
#include "n_piPlus_from_root.h"
//#include "MM_root_analysis.h"
#include "MM_RootDataAnalysis.h"

int main() {

    // ChoiceRoot test("pass2.txt", "output.root");
    // test.analysisCycle();

    // Neutron_PiPlus_Root_4D_MM test("output.root");
    // test.analysisCycle();
    // test.fitingResults();

    //MM_Root_Analysis test(36);
    // test.handFitting(0, 1);
    //test.prepareResults_25_FullFit();

    // MM_FullFitAnalysis test;
    // test.analysisCycle();

    // MM_ShortFitAnalysis_4D test1("fitShortMM_25.root");
    // test1.analysisCycle();

    // MM_FittedShortDataAnalysis_4D test2("fitShortMM_25.root");
    // test2.analysisCycle();

    // Neutron_PiPlus_Root_4D_MM test1("output.root");
    // test1.analysisCycle();

    // MM_FullFitAnalysis_4D test2("fitMM_50.root");
    // test2.analysisCycle();

    // MM_FittedDataAnalysis_4D test3("fitMM_50.root");
    // test3.analysisCycle();

    // MM_FittedDataAnalysis_4D_PRESENTATION
    // MM_FittedDataAnalysis_4D_PRESENTATION test3("fitMM_50.root");
    // test3.analysisCycle();






    // SortedOutputCreator test0("output.root");
    // test0.analysisCycle();

    // MM_FullSortedDataAnalysis_4D test1;
    // test1.analysisCycle();

    // MM_FittingSortedDataAnalysis_4D test2("fitMM_var50.root");
    // test2.analysisCycle();

    // MM_DrawFittedShortSortedDataAnalysis_4D test3("fitMM_var50.root");
    // test3.analysisCycle();

    MM_FittedSortedDataAnalysis_4D_PRESENTATION test4("fitMM_var50.root");
    test4.analysisCycle();

    return 0;
}