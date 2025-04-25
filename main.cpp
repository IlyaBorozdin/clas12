#include <string>
#include <iostream>
#include <fstream>
#include <vector>

//#include "neutron_piPlus.h"
//#include "yield_in_q2-w.h"
//#include "choiceRoot.h"
// #include "n_piPlus_from_root.h"
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

    MM_FullFitAnalysis_4D test;
    test.analysisCycle();

    return 0;
}