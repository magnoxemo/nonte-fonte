#ifndef NONTE_FONTE_FET_ENGINE_H
#define NONTE_FONTE_FET_ENGINE_H

#include "Legendre.h"

#include <omp.h>
#include <random>
#include <stdexcept>
#include <vector>

namespace nontefonte {

class Function;

/**
 * @brief Base class for the Fonte Monte Carlo functional expansion engine.
 *
 * Provides a dimension-agnostic framework for estimating probability density
 * functions via Monte Carlo sampling and Legendre polynomial expansions.
 * Concrete subclasses (1-D, 2-D, 3-D) fix the dimensionality and must
 * implement @ref evaluatePDF and @ref basisProduct, and may override
 * @ref accumulateTrials for dimension-specific optimisations.
 *
 * @note All public data members are intentionally exposed for performance;
 *       prefer @ref getCoefficients() for read access to the result.
 */
class FonteEngineBase {
public:
  /**
   * @brief Constructs the engine with a default domain of [-1, 1]^N.
   *
   * @param orders           Polynomial expansion order for each dimension.
   *                         The length of this vector defines the
   * dimensionality.
   * @param number_of_trials Number of Monte Carlo samples to draw per
   * simulation.
   */
  FonteEngineBase(std::vector<int> orders, long int number_of_trials);

  /**
   * @brief Constructs the engine with an explicit integration domain.
   *
   * @param orders           Polynomial expansion order for each dimension.
   * @param number_of_trials Number of Monte Carlo samples to draw per
   * simulation.
   * @param domain           Integration bounds for each dimension as
   *                         {min, max} pairs. Must have the same length as @p
   * orders.
   */
  FonteEngineBase(std::vector<int> orders, long int number_of_trials,
                  std::vector<std::pair<double, double>> domain);

  /** @brief Virtual destructor for safe polymorphic deletion. */
  virtual ~FonteEngineBase() = default;

  /**
   * @brief Runs the Monte Carlo simulation and populates @ref _co_efficients.
   *
   * Spawns one OpenMP thread per available worker. Each thread draws an equal
   * share of @ref _number_of_trials sample points, evaluates the PDF and
   * Legendre basis products, and accumulates partial sums that are reduced
   * into @ref _co_efficients at the end.
   *
   * @param pdf The probability density function to approximate.
   */
  void runSimulation(Function &pdf);

  /**
   * @brief Converts a multi-dimensional index to a flat (linear) array index.
   *
   * Uses a row-major (C-order) stride calculation so that the last dimension
   * varies fastest.
   *
   * @param multi_index Per-dimension indices. Must have the same length as
   *                    @ref _orders.
   * @return The corresponding flat index into @ref _co_efficients.
   */
  int flattenIndex(const std::vector<int> &multi_index) const {
    int index = 0, stride = 1;
    for (int d = static_cast<int>(_orders.size()) - 1; d >= 0; --d) {
      index += multi_index[d] * stride;
      stride *= (_orders[d] + 1);
    }
    return index;
  }

  /**
   * @brief Returns a copy of the computed expansion coefficients.
   *
   * Call this after @ref runSimulation() to retrieve the Legendre
   * expansion coefficients c_{i,j,...} that approximate the PDF.
   *
   * @return A vector of length @ref _number_of_co_efficients containing
   *         the expansion coefficients in row-major order.
   */
  std::vector<double> getCoefficients() const { return _co_efficients; }

  // --- data ---

  /** @brief Total number of Monte Carlo trials per simulation run. */
  const long int _number_of_trials;

  /** @brief Total number of expansion coefficients (product of all (order+1)).
   */
  int _number_of_co_efficients;

  /** @brief Polynomial expansion order for each dimension. */
  std::vector<int> _orders;

  /** @brief Computed Legendre expansion coefficients; populated by @ref
   * runSimulation(). */
  std::vector<double> _co_efficients;

  /** @brief Integration domain as {min, max} pairs, one per dimension. */
  std::vector<std::pair<double, double>> _domain;

protected:
  /**
   * @brief Evaluates the PDF at a given sample point.
   *
   * Pure virtual — each subclass must implement this to forward the call
   * to @p pdf with the correct number of arguments for its dimensionality.
   *
   * @param pdf The probability density function.
   * @param x   Sample point in the physical domain, one value per dimension.
   * @return    The PDF value at @p x.
   */
  virtual double evaluatePDF(Function &pdf,
                             const std::vector<double> &x) const = 0;

  /**
   * @brief Computes the product of Legendre basis polynomials at a sample
   * point.
   *
   * Pure virtual — each subclass must implement this to compute:
   * @f[ \prod_{d} P_{mi_d}(x_d) @f]
   * where @f$ P_n @f$ is the n-th Legendre polynomial. Subclasses with a
   * fixed dimensionality can unroll the product loop for performance.
   *
   * @param basis Per-dimension Legendre basis objects.
   * @param x     Sample point in the physical domain.
   * @param mi    Multi-index specifying the polynomial order per dimension.
   * @return      The scalar basis product at @p x.
   */
  virtual double basisProduct(const std::vector<LegendreBasis> &basis,
                              const std::vector<double> &x,
                              const std::vector<int> &mi) const = 0;

  /**
   * @brief Runs all Monte Carlo trials and returns the accumulated coefficient
   * sums.
   *
   * The base implementation runs a sequential loop over @p trials_this_thread
   * sample points, calling @ref evaluatePDF and @ref basisProduct (both
   * provided by the subclass) at each point. Subclasses may override this
   * to apply further optimisations, but typically the two pure virtuals
   * are sufficient to specialise behaviour.
   *
   * @param pdf              The probability density function (thread-local
   * copy).
   * @param basis            Per-dimension Legendre basis objects.
   * @param dists            Per-dimension uniform distributions over the
   * domain.
   * @param gen              Thread-local random number engine.
   * @param trials_this_thread Number of trials this thread is responsible for.
   * @return A vector of length @ref _number_of_co_efficients containing the
   *         un-normalised partial coefficient sums for this thread.
   */
  virtual std::vector<double>
  accumulateTrials(Function &pdf, std::vector<LegendreBasis> &basis,
                   std::vector<std::uniform_real_distribution<double>> &dists,
                   std::default_random_engine &gen,
                   long trials_this_thread) const;

  /**
   * @brief Returns the number of dimensions this engine expects.
   *
   * Implemented by each concrete subclass to return 1, 2, or 3.
   * Used by @ref checkDimensions() to validate constructor arguments.
   *
   * @return The required dimensionality of this engine.
   */
  virtual int expectedDimensions() const = 0;

  /**
   * @brief Validates that @ref _orders has the correct number of dimensions.
   *
   * Should be called at the end of each concrete subclass constructor.
   *
   * @throws std::invalid_argument if the length of @ref _orders does not
   *         equal @ref expectedDimensions().
   */
  void checkDimensions() const {
    if (static_cast<int>(_orders.size()) != expectedDimensions())
      throw std::invalid_argument(
          "Order vector length does not match engine dimensionality.");
  }

  /**
   * @brief Pre-computed table of multi-indices, indexed by flat coefficient
   * index.
   *
   * Built once in @ref initCoefficients() so that @ref accumulateTrials
   * never allocates a temporary index vector inside the hot loop.
   */
  std::vector<std::vector<int>> _multi_index_table;

private:
  /**
   * @brief Allocates and zero-initialises @ref _co_efficients and
   *        populates @ref _multi_index_table.
   *
   * Computes @ref _number_of_co_efficients as the product of all
   * (@ref _orders[d] + 1), resizes the coefficient vector, and
   * pre-computes every multi-index entry.
   */
  void initCoefficients();

  /**
   * @brief Converts a flat index back into a multi-dimensional index.
   *
   * The inverse of @ref flattenIndex(). Used only during initialisation
   * to populate @ref _multi_index_table.
   *
   * @param flat A linear index in [0, @ref _number_of_co_efficients).
   * @return     The corresponding per-dimension index vector.
   */
  std::vector<int> multiIndex(int flat) const;

  /**
   * @brief Draws a single sample point from the integration domain.
   *
   * Samples one value from each per-dimension distribution and returns
   * them as a coordinate vector.
   *
   * @param dists Per-dimension uniform distributions.
   * @param gen   Random number engine.
   * @return      A sample point with one coordinate per dimension.
   */
  std::vector<double>
  samplePoint(std::vector<std::uniform_real_distribution<double>> &dists,
              std::default_random_engine &gen) const;
};

} // namespace nontefonte

#endif