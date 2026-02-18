#include <random>
#include <omp.h>
#include <iostream>

#include "FETEngine.h"
#include "PDF.h"



int main(int argc, char* argv[]){

    const long int N = std::stoi(argv[1]);
    float average = 0.0f;
    nontefonte::Function pdf("sin(x)/cos(y)+(y-x/cos(x))^2 - sqrt(x*x)", {"x", "y"});

    nontefonte::FETEngine fet_engine({8,8}, N);
    fet_engine.runSimulation(pdf);

    for (int i = 0; i<fet_engine._orders[0]; i++){
        for (int j = 0; j<fet_engine._orders[1]; j++){
            std::cout<<fet_engine._co_efficients[i* fet_engine._orders[0] +j]<<" ";
        }
        printf("\n");
    }

    return 0;

}
