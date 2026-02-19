#ifndef NONTE_FONTE_FONTE_ENGINE_TWO_DIMENSION_H
#define NONTE_FONTE_FONTE_ENGINE_TWO_DIMENSION_H

#include "FonteEngineBase.h"

namespace nontefonte {

    /**
     * @brief 2-D specialisation of the Fonte Monte Carlo expansion engine.
     *
     * Fixes the dimensionality to 2 and provides optimised overrides of the
     * pure virtual methods from @ref FonteEngineBase, unrolling all dimension
     * loops to direct two-term expressions.
     */
    class FonteEngineTwoDimension : public FonteEngineBase {
    public:

        /**
         * @brief Constructs a 2-D engine with the default domain [-1, 1]^2.
         *
         * @param orders           Polynomial expansion orders as {order_x, order_y}.
         * @param number_of_trials Number of Monte Carlo samples per simulation.
         */
        FonteEngineTwoDimension(std::pair<int, int> orders,
                                long int number_of_trials);

        /**
         * @brief Constructs a 2-D engine with an explicit integration domain.
         *
         * @param orders           Polynomial expansion orders as {order_x, order_y}.
         * @param number_of_trials Number of Monte Carlo samples per simulation.
         * @param domain           Integration bounds for each dimension as
         *                         {min, max} pairs. Must have exactly 2 elements.
         */
        FonteEngineTwoDimension(std::pair<int, int> orders,
                                long int number_of_trials,
                                std::vector<std::pair<double, double>> domain);

    protected:

        /**
         * @brief Returns 2, the fixed dimensionality of this engine.
         * @return 2
         */
        int expectedDimensions() const override;

        /**
         * @brief Evaluates the 2-D PDF by calling @p pdf with two arguments.
         *
         * @param pdf The probability density function.
         * @param x   Sample point; `x[0]` and `x[1]` are used.
         * @return    The PDF value at `(x[0], x[1])`.
         */
        double evaluatePDF(Function &pdf,
                           const std::vector<double> &x) const override;

        /**
         * @brief Computes the 2-D Legendre basis product at a sample point.
         *
         * Returns the direct two-term product:
         * @f[ P_{mi_0}(x_0) \cdot P_{mi_1}(x_1) @f]
         *
         * @param basis Two-element basis vector.
         * @param x     Sample point; `x[0]` and `x[1]` are used.
         * @param mi    Multi-index; `mi[0]` and `mi[1]` are used.
         * @return      The scalar basis product at `(x[0], x[1])`.
         */
        double basisProduct(const std::vector<LegendreBasis> &basis,
                            const std::vector<double> &x,
                            const std::vector<int> &mi) const override;
    };

} // namespace nontefonte

#endif // NONTE_FONTE_FONTE_ENGINE_TWO_DIMENSION_H