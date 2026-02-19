#ifndef NONTE_FONTE_FONTE_ENGINE_ONE_DIMENSION_H
#define NONTE_FONTE_FONTE_ENGINE_ONE_DIMENSION_H

#include "FonteEngineBase.h"

namespace nontefonte {

/**
 * @brief 1-D specialisation of the Fonte Monte Carlo expansion engine.
 *
 * Fixes the dimensionality to 1 and provides optimised overrides of the
 * pure virtual methods from @ref FonteEngineBase that avoid any inner
 * loop or runtime dimension check.
 */
class FonteEngineOneDimension : public FonteEngineBase {
public:
  /**
   * @brief Constructs a 1-D engine with the default domain [-1, 1].
   *
   * @param order            Polynomial expansion order.
   * @param number_of_trials Number of Monte Carlo samples per simulation.
   */
  FonteEngineOneDimension(int order, long int number_of_trials);

  /**
   * @brief Constructs a 1-D engine with an explicit integration domain.
   *
   * @param order            Polynomial expansion order.
   * @param number_of_trials Number of Monte Carlo samples per simulation.
   * @param domain           Integration bounds as a {min, max} pair.
   */
  FonteEngineOneDimension(int order, long int number_of_trials,
                          std::pair<double, double> domain);

protected:
  /**
   * @brief Returns 1, the fixed dimensionality of this engine.
   * @return 1
   */
  int expectedDimensions() const override;

  /**
   * @brief Evaluates the 1-D PDF by calling @p pdf with a single argument.
   *
   * @param pdf The probability density function.
   * @param x   Sample point; only `x[0]` is used.
   * @return    The PDF value at `x[0]`.
   */
  double evaluatePDF(Function &pdf,
                     const std::vector<double> &x) const override;

  /**
   * @brief Computes the 1-D Legendre basis value at a sample point.
   *
   * @param basis Single-element basis vector; only `basis[0]` is used.
   * @param x     Sample point; only `x[0]` is used.
   * @param mi    Multi-index; only `mi[0]` is used.
   * @return      @f$ P_{mi_0}(x_0) @f$
   */
  double basisProduct(const std::vector<LegendreBasis> &basis,
                      const std::vector<double> &x,
                      const std::vector<int> &mi) const override;
};

} // namespace nontefonte

#endif // NONTE_FONTE_FONTE_ENGINE_ONE_DIMENSION_H