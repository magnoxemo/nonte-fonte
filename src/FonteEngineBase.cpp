#include "FonteEngineBase.h"
#include "PDF.h"

#include <algorithm>
#include <omp.h>
#include <random>
#include <stdexcept>

namespace nontefonte {

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

    FonteEngineBase::FonteEngineBase(std::vector<int> orders,
                                     long int number_of_trials)
            : _number_of_trials(number_of_trials), _orders(std::move(orders)) {
      _domain.reserve(_orders.size());
      for (size_t i = 0; i < _orders.size(); ++i)
        _domain.emplace_back(-1.0, 1.0);
      initCoefficients();
    }

    FonteEngineBase::FonteEngineBase(std::vector<int> orders,
                                     long int number_of_trials,
                                     std::vector<std::pair<double, double>> domain)
            : _number_of_trials(number_of_trials),
              _orders(std::move(orders)),
              _domain(std::move(domain)) {
      initCoefficients();
    }

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

    void FonteEngineBase::initCoefficients() {
      int total = 1;
      for (int o : _orders)
        total *= (o + 1);

      _number_of_co_efficients = total;
      _co_efficients.assign(total, 0.0);

      // Pre-compute multi-index and normalization tables once at construction.
      // This keeps the hot trial loop free of any per-iteration allocation or
      // repeated arithmetic.
      _multi_index_table.resize(total);
      _norm_factors.resize(total);

      for (int flat = 0; flat < total; ++flat) {
        _multi_index_table[flat] = multiIndex(flat);

        const auto &mi = _multi_index_table[flat];
        double factor = 1.0;
        for (size_t d = 0; d < mi.size(); ++d) {
          const double width = _domain[d].second - _domain[d].first;
          /* Just reexpression my thoughts:
           * For montecarlo simulation we use this formula as the normalization coefficient
           * (2 * n + 1) /2 which is formulated presuming rho(x) = 1
           * But when we are sampling random PDF the weighting function is basically rho = 1 /(b -a)
           * random number related probability so actual normalization coefficient will be
           * (2 * n + 1) */
          factor *= static_cast<double>(2 * mi[d] + 1) ;
        }
        _norm_factors[flat] = factor;
      }
    }

    std::vector<int> FonteEngineBase::multiIndex(int flat) const {
      const int dim = static_cast<int>(_orders.size());
      std::vector<int> idx(dim);
      for (int d = dim - 1; d >= 0; --d) {
        idx[d] = flat % (_orders[d] + 1);
        flat   /= (_orders[d] + 1);
      }
      return idx;
    }

    std::vector<double> FonteEngineBase::samplePoint(
            std::vector<std::uniform_real_distribution<double>> &dists,
            std::default_random_engine &gen) const {
      std::vector<double> x(_orders.size());
      for (size_t d = 0; d < x.size(); ++d)
        x[d] = dists[d](gen);
      return x;
    }

// ---------------------------------------------------------------------------
// accumulateTrials
//
// Runs this thread's share of Monte Carlo trials and returns un-normalised
// partial sums. The (2n+1)/(b-a) normalization is applied here per
// coefficient before returning, so runSimulation only needs to divide by N
// and reduce into the shared array.
// ---------------------------------------------------------------------------

    std::vector<double> FonteEngineBase::accumulateTrials(
            Function &pdf,
            std::vector<LegendreBasis> &basis,
            std::vector<std::uniform_real_distribution<double>> &dists,
            std::default_random_engine &gen,
            long trials_this_thread) const {

      std::vector<double> partial_sums(_number_of_co_efficients, 0.0);

      for (long i = 0; i < trials_this_thread; ++i) {
        const auto   x       = samplePoint(dists, gen);
        const double pdf_val = evaluatePDF(pdf, x);

        for (int flat = 0; flat < _number_of_co_efficients; ++flat) {
          partial_sums[flat] +=
                  basisProduct(basis, x, _multi_index_table[flat]) * pdf_val;
        }
      }

      // Apply the Legendre normalization factor per coefficient.
      // The remaining 1/N scaling is applied in runSimulation after reduction.
      for (int flat = 0; flat < _number_of_co_efficients; ++flat)
        partial_sums[flat] *= _norm_factors[flat];

      return partial_sums;
    }

// ---------------------------------------------------------------------------
// runSimulation
// ---------------------------------------------------------------------------

    void FonteEngineBase::runSimulation(Function &pdf) {
      const int dim = static_cast<int>(_orders.size());

      std::fill(_co_efficients.begin(), _co_efficients.end(), 0.0);

#pragma omp parallel
      {
        const int  num_threads        = omp_get_num_threads();
        const int  thread_id          = omp_get_thread_num();
        const long base_trials        = _number_of_trials / num_threads;
        const long remainder          = _number_of_trials % num_threads;
        const long trials_this_thread =
                base_trials + (thread_id == 0 ? remainder : 0L);

        std::default_random_engine gen(
                std::random_device{}() + static_cast<unsigned>(thread_id));

        Function local_pdf = pdf;

        std::vector<LegendreBasis> basis;
        std::vector<std::uniform_real_distribution<double>> dists;
        basis.reserve(dim);
        dists.reserve(dim);
        for (int d = 0; d < dim; ++d) {
          basis.emplace_back(_orders[d]);
          dists.emplace_back(_domain[d].first, _domain[d].second);
        }

        auto partial_sums =
                accumulateTrials(local_pdf, basis, dists, gen, trials_this_thread);

        // Divide by total N — normalization factor is already baked in.
        const double inv_total = 1.0 / static_cast<double>(_number_of_trials);
        for (double &s : partial_sums)
          s *= inv_total;

#pragma omp critical
        for (int flat = 0; flat < _number_of_co_efficients; ++flat)
          _co_efficients[flat] += partial_sums[flat];
      }
    }

} // namespace nontefonte