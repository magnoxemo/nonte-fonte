#include <iostream>

#include "FonteEngineOneDimension.h"
#include "PDF.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <number_of_trials>\n";
    return 1;
  }

  const long int number_of_trials = std::stol(argv[1]);

  nontefonte::Function pdf("sin(x)+cos(x)^2", {"x"});

  nontefonte::FonteEngineOneDimension engine(8, number_of_trials);
  engine.runSimulation(pdf);

  const auto coefficients = engine.getCoefficients();
  const int order_x = engine._orders[0];

  std::cout << "[ ";
  for (int i = 0; i <= order_x; ++i) {
      std::cout << coefficients[i ] << " ,";
  }
  std::cout << "]\n";

  return 0;
}