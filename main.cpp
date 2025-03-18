#include "neutron_piPlus.h"
//#include "choiceRoot.h"

int main() {
    
    vector<string> fileNames = {
        "rec_clas_005893.evio.00315-00319.hipo",
        "rec_clas_005893.evio.00390-00394.hipo",
        "rec_clas_005893.evio.00425-00429.hipo",
        "rec_clas_005893.evio.00470-00474.hipo",
    };
    
    
    Neutron_PiPlus_Detailed_Analysis_MM test(fileNames);
    //ChoiceRoot test("input.txt", "output.root");

    test.analysisCycle();

    return 0;
}