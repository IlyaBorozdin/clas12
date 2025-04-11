#include <string>
#include <iostream>
#include <fstream>
#include <vector>

//#include "neutron_piPlus.h"
//#include "yield_in_q2-w.h"
//#include "choiceRoot.h"
// #include "n_piPlus_from_root.h"
#include "MM_root_analysis.h"

int main() {

    /*
    std::string filename = "pass2.txt";
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл: " << filename << std::endl;
        return 1;
    }

    std::vector<std::string> filenames;
    std::string path;

    while (std::getline(file, path)) {
        filenames.push_back(path);
    }

    file.close();

    std::cout << "Прочитанные пути из файла:" << std::endl;
    for (const auto& path : filenames) {
        std::cout << path << std::endl;
    }
    
    Neutron_PiPlus_Yield_Analysis_MM test(filenames);

    test.analysisCycle();
    */

    // ChoiceRoot test("pass2.txt", "output.root");
    // test.analysisCycle();

    // Neutron_PiPlus_Root_MM test("output.root");
    // test.analysisCycle();
    // test.fitingResults();

    MM_Root_Analysis test(36);
    // test.handFitting(0, 1);
    test.prepareResults_25_FullFit();

    return 0;
}