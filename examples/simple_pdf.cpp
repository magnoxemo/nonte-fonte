#include <iostream>

#include "FonteEngineThreeDimension.h"
#include "PDF.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <number_of_trials>\n";
    return 1;
  }

  const long int number_of_trials = std::stol(argv[1]);

  nontefonte::Function pdf("sin(x)/cos(y)+(y-x/cos(x))^2 +z*z - sqrt(x*x)",
                           {"x", "y", "z"});

  nontefonte::FonteEngineThreeDimension engine({8, 8, 9}, number_of_trials);
  engine.runSimulation(pdf);

  const auto coefficients = engine.getCoefficients();
  const int order_x = engine._orders[0];
  const int order_y = engine._orders[1];

  for (int i = 0; i <= order_x; ++i) {
    for (int j = 0; j <= order_y; ++j) {
      std::cout << coefficients[i * (order_y + 1) + j] << " ";
    }
    std::cout << "\n";
  }

  return 0;
}