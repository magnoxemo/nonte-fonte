#include <random>
#include <omp.h>
#include <iostream>

#include "FETEngine.h"
#include "PDF.h"

int main(int argc, char* argv[]) {
    const long int N = std::stoi(argv[1]);

    nontefonte::Function pdf("sin(x)/cos(y)+(y-x/cos(x))^2 +z*z - sqrt(x*x)", {"x", "y", "z"});

    nontefonte::FETEngine fet_engine({8, 8, 9}, N);
    fet_engine.runSimulation(pdf);

    for (int i = 0; i <= fet_engine._orders[0]; i++) {
        for (int j = 0; j <= fet_engine._orders[1]; j++) {
            std::cout << fet_engine._co_efficients[i * (fet_engine._orders[1] + 1) + j] << " ";
        }
        std::cout << "\n";
    }

    return 0;
}