#include "FonteEngineOneDimension.h"
#include "PDF.h"

namespace nontefonte {

FonteEngineOneDimension::FonteEngineOneDimension(int order,
                                                 long int number_of_trials)
    : FonteEngineBase({order}, number_of_trials) {
  checkDimensions();
}

FonteEngineOneDimension::FonteEngineOneDimension(
    int order, long int number_of_trials, std::pair<double, double> domain)
    : FonteEngineBase({order}, number_of_trials, {domain}) {
  checkDimensions();
}

int FonteEngineOneDimension::expectedDimensions() const { return 1; }

double
FonteEngineOneDimension::evaluatePDF(Function &pdf,
                                     const std::vector<double> &x) const {
  return pdf(x[0]);
}

double
FonteEngineOneDimension::basisProduct(const std::vector<LegendreBasis> &basis,
                                      const std::vector<double> &x,
                                      const std::vector<int> &mi) const {
  return basis[0](x[0], mi[0]);
}

} // namespace nontefonte