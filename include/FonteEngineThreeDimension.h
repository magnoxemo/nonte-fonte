#ifndef NONTE_FONTE_FONTE_ENGINE_THREE_DIMENSION_H
#define NONTE_FONTE_FONTE_ENGINE_THREE_DIMENSION_H

#include "FonteEngineBase.h"

#include <tuple>

namespace nontefonte {

    /**
     * @brief 3-D specialisation of the Fonte Monte Carlo expansion engine.
     *
     * Fixes the dimensionality to 3 and provides optimised overrides of the
     * pure virtual methods from @ref FonteEngineBase, unrolling all dimension
     * loops to direct three-term expressions.
     */
    class FonteEngineThreeDimension : public FonteEngineBase {
    public:

        /**
         * @brief Constructs a 3-D engine with the default domain [-1, 1]^3.
         *
         * @param orders           Polynomial expansion orders as
         *                         {order_x, order_y, order_z}.
         * @param number_of_trials Number of Monte Carlo samples per simulation.
         */
        FonteEngineThreeDimension(std::tuple<int, int, int> orders,
                                  long int number_of_trials);

        /**
         * @brief Constructs a 3-D engine with an explicit integration domain.
         *
         * @param orders           Polynomial expansion orders as
         *                         {order_x, order_y, order_z}.
         * @param number_of_trials Number of Monte Carlo samples per simulation.
         * @param domain           Integration bounds for each dimension as
         *                         {min, max} pairs. Must have exactly 3 elements.
         */
        FonteEngineThreeDimension(std::tuple<int, int, int> orders,
                                  long int number_of_trials,
                                  std::vector<std::pair<double, double>> domain);

    protected:

        /**
         * @brief Returns 3, the fixed dimensionality of this engine.
         * @return 3
         */
        int expectedDimensions() const override;

        /**
         * @brief Evaluates the 3-D PDF by calling @p pdf with three arguments.
         *
         * @param pdf The probability density function.
         * @param x   Sample point; `x[0]`, `x[1]`, and `x[2]` are used.
         * @return    The PDF value at `(x[0], x[1], x[2])`.
         */
        double evaluatePDF(Function &pdf,
                           const std::vector<double> &x) const override;

        /**
         * @brief Computes the 3-D Legendre basis product at a sample point.
         *
         * Returns the direct three-term product:
         * @f[ P_{mi_0}(x_0) \cdot P_{mi_1}(x_1) \cdot P_{mi_2}(x_2) @f]
         *
         * @param basis Three-element basis vector.
         * @param x     Sample point; `x[0]`, `x[1]`, and `x[2]` are used.
         * @param mi    Multi-index; `mi[0]`, `mi[1]`, and `mi[2]` are used.
         * @return      The scalar basis product at `(x[0], x[1], x[2])`.
         */
        double basisProduct(const std::vector<LegendreBasis> &basis,
                            const std::vector<double> &x,
                            const std::vector<int> &mi) const override;
    };

} // namespace nontefonte

#endif // NONTE_FONTE_FONTE_ENGINE_THREE_DIMENSION_H