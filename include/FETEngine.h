#ifndef NONTE_FONTE_FET_ENGINE_H
#define NONTE_FONTE_FET_ENGINE_H

#include "Legendre.h"
#include <omp.h>
#include <random>
#include <vector>

namespace nontefonte {



// forward declaration
class Function;

class FETEngine {

public:
  FETEngine(std::vector<int> orders, long int number_of_trials);
  FETEngine(std::vector<int> orders, long int number_of_trials,
            std::vector<std::pair<double, double>> domain);

  void runSimulation(Function &pdf);

  int flattenIndex(const std::vector<int> &multi_index) const {
    int index = 0;
    int stride = 1;
    for (int d = _orders.size() - 1; d >= 0; --d) {
      index += multi_index[d] * stride;
      stride *= (_orders[d] + 1);
    }
    return index;
  }

  const long int _number_of_trials;
  int _number_of_co_efficients;

  std::vector<std::pair<double, double>> _domain;
  std::vector<int> _orders;
  std::vector<double> _co_efficients;
};
} // namespace nontefonte

#endif // NONTE_FONTE_FETENGINE_H
