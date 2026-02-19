#include "FonteEngineThreeDimension.h"
#include "PDF.h"

namespace nontefonte {

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

    FonteEngineThreeDimension::FonteEngineThreeDimension(
            std::tuple<int, int, int> orders,
            long int number_of_trials)
            : FonteEngineBase(
            {std::get<0>(orders), std::get<1>(orders), std::get<2>(orders)},
            number_of_trials) {
        checkDimensions();
    }

    FonteEngineThreeDimension::FonteEngineThreeDimension(
            std::tuple<int, int, int> orders,
            long int number_of_trials,
            std::vector<std::pair<double, double>> domain)
            : FonteEngineBase(
            {std::get<0>(orders), std::get<1>(orders), std::get<2>(orders)},
            number_of_trials, std::move(domain)) {
        checkDimensions();
    }

// ---------------------------------------------------------------------------
// Protected overrides
// ---------------------------------------------------------------------------

    int FonteEngineThreeDimension::expectedDimensions() const {
        return 3;
    }

    double FonteEngineThreeDimension::evaluatePDF(Function &pdf,
                                                  const std::vector<double> &x) const {
        return pdf(x[0], x[1], x[2]);
    }

    double FonteEngineThreeDimension::basisProduct(
            const std::vector<LegendreBasis> &basis,
            const std::vector<double> &x,
            const std::vector<int> &mi) const {
        return basis[0](x[0], mi[0])
               * basis[1](x[1], mi[1])
               * basis[2](x[2], mi[2]);
    }

} // namespace nontefonte