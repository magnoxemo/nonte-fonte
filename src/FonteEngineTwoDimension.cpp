#include "FonteEngineTwoDimension.h"
#include "PDF.h"

namespace nontefonte {

FonteEngineTwoDimension::FonteEngineTwoDimension(std::pair<int, int> orders,
                                                 long int number_of_trials)
    : FonteEngineBase({orders.first, orders.second}, number_of_trials) {
  checkDimensions();
}

FonteEngineTwoDimension::FonteEngineTwoDimension(
    std::pair<int, int> orders, long int number_of_trials,
    std::vector<std::pair<double, double>> domain)
    : FonteEngineBase({orders.first, orders.second}, number_of_trials,
                      std::move(domain)) {
  checkDimensions();
}


int FonteEngineTwoDimension::expectedDimensions() const { return 2; }

double
FonteEngineTwoDimension::evaluatePDF(Function &pdf,
                                     const std::vector<double> &x) const {
  return pdf(x[0], x[1]);
}

double
FonteEngineTwoDimension::basisProduct(const std::vector<LegendreBasis> &basis,
                                      const std::vector<double> &x,
                                      const std::vector<int> &mi) const {
  return basis[0](x[0], mi[0]) * basis[1](x[1], mi[1]);
}

} // namespace nontefonte