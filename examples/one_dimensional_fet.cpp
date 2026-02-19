#include <iostream>
#include <vector>

#include "FonteEngineOneDimension.h"
#include "PDF.h"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <number_of_trials> <order>\n";
    return 1;
  }

  const long number_of_trials = std::stol(argv[1]);
  const int order = std::stoi(argv[2]);

  nontefonte::Function pdf("3*x*x-4*x+5", {"x"});

  nontefonte::FonteEngineOneDimension engine(order, number_of_trials);
  engine.runSimulation(pdf);

  const auto coefficients = engine.getCoefficients();

  std::cout << "[";
  for (int i = 0; i <= order; ++i) {
    std::cout << coefficients[i];
    if (i < order)
      std::cout << ", ";
  }
  std::cout << "]\n";

  return 0;
}