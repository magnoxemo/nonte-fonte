#ifndef NONTE_FONTE_FET_ENGINE_H
#define NONTE_FONTE_FET_ENGINE_H

#include "Legendre.h"

#include <omp.h>
#include <random>
#include <vector>

namespace nontefonte {

class Function;

class FETEngine {
public:
  FETEngine(std::vector<int> orders, long int number_of_trials);
  FETEngine(std::vector<int> orders, long int number_of_trials,
            std::vector<std::pair<double, double>> domain);

  void runSimulation(Function &pdf);

  int flattenIndex(const std::vector<int> &multi_index) const {
    int index = 0, stride = 1;
    for (int d = _orders.size() - 1; d >= 0; --d) {
      index += multi_index[d] * stride;
      stride *= (_orders[d] + 1);
    }
    return index;
  }

  // public data (consider getters if this grows)
  const long int _number_of_trials;
  int _number_of_co_efficients;
  std::vector<int> _orders;
  std::vector<double> _co_efficients;
  std::vector<std::pair<double, double>> _domain;

  std::vector<double> getCoefficients(){
        return _co_efficients;
  };

private:
  void initCoefficients();
  std::vector<int> multiIndex(int flat) const;

  std::vector<double>
  samplePoint(std::vector<std::uniform_real_distribution<double>> &dists,
              std::default_random_engine &gen) const;

  double evaluatePDF(Function &pdf, const std::vector<double> &x) const;

  double basisProduct(const std::vector<LegendreBasis> &basis,
                      const std::vector<double> &x,
                      const std::vector<int> &mi) const;

  std::vector<double>
  accumulateTrials(Function &pdf, std::vector<LegendreBasis> &basis,
                   std::vector<std::uniform_real_distribution<double>> &dists,
                   std::default_random_engine &gen) const;
};

} // namespace nontefonte

#endif