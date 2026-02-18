#include <iostream>
#include <vector>

#include "FETEngine.h"
#include "PDF.h"

int main(int argc, char* argv[]) {
    const long N    = std::stol(argv[1]);
    const int  order = std::stoi(argv[2]);

    nontefonte::Function pdf("sin(15*x)", {"x"});

    nontefonte::FETEngine engine({order}, N);
    engine.runSimulation(pdf);

    const auto& coefficients = engine.getCoefficients();

    std::cout << "[";
    for (int i = 0; i <= order; i++) {
        std::cout << coefficients[i];
        if (i < order) std::cout << ", ";
    }
    std::cout << "]\n";
}