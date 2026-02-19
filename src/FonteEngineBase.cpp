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
    : _number_of_trials(number_of_trials), _orders(std::move(orders)),
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

  // Pre-compute every multi-index once so the hot trial loop can index
  // directly into this table without any per-iteration allocation.
  _multi_index_table.resize(total);
  for (int flat = 0; flat < total; ++flat)
    _multi_index_table[flat] = multiIndex(flat);
}

std::vector<int> FonteEngineBase::multiIndex(int flat) const {
  const int dim = static_cast<int>(_orders.size());
  std::vector<int> idx(dim);
  for (int d = dim - 1; d >= 0; --d) {
    idx[d] = flat % (_orders[d] + 1);
    flat /= (_orders[d] + 1);
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
// accumulateTrials  (non-pure default — subclasses may override)
//
// Runs a self-contained sequential loop over this thread's share of trials.
// evaluatePDF and basisProduct are pure virtual, so the subclass controls
// how the PDF is called and how the basis product is computed; this method
// only owns the loop structure and accumulation logic.
// ---------------------------------------------------------------------------

std::vector<double> FonteEngineBase::accumulateTrials(
    Function &pdf, std::vector<LegendreBasis> &basis,
    std::vector<std::uniform_real_distribution<double>> &dists,
    std::default_random_engine &gen, long trials_this_thread) const {

  std::vector<double> partial_sums(_number_of_co_efficients, 0.0);

  for (long i = 0; i < trials_this_thread; ++i) {
    const auto x = samplePoint(dists, gen);
    const double pdf_val = evaluatePDF(pdf, x); // subclass-provided

    for (int flat = 0; flat < _number_of_co_efficients; ++flat) {
      // _multi_index_table avoids a vector allocation per iteration.
      partial_sums[flat] +=
          basisProduct(basis, x, _multi_index_table[flat]) * pdf_val;
    }
  }
  return partial_sums;
}

// ---------------------------------------------------------------------------
// runSimulation
//
// Owns the parallel decomposition. Each thread:
//   1. Builds its own RNG, PDF copy, basis evaluators and distributions.
//   2. Calls accumulateTrials for its share of the total trial count.
//   3. Normalizes its partial sums before entering the critical section,
//      keeping the shared reduction as cheap as possible.
// ---------------------------------------------------------------------------

void FonteEngineBase::runSimulation(Function &pdf) {
  const int dim = static_cast<int>(_orders.size());

  // Zero coefficients so runSimulation() is safe to call more than once.
  std::fill(_co_efficients.begin(), _co_efficients.end(), 0.0);

#pragma omp parallel
  {
    const int num_threads = omp_get_num_threads();
    const int thread_id = omp_get_thread_num();

    // Distribute trials as evenly as possible; thread 0 absorbs the
    // remainder so the total is always exactly _number_of_trials.
    const long base_trials = _number_of_trials / num_threads;
    const long remainder = _number_of_trials % num_threads;
    const long trials_this_thread =
        base_trials + (thread_id == 0 ? remainder : 0L);

    // Unique per-thread seed: blends entropy from random_device with
    // the thread index to guarantee distinct RNG streams.
    std::default_random_engine gen(std::random_device{}() +
                                   static_cast<unsigned>(thread_id));

    // Thread-local copy of the PDF functor — required when the functor
    // carries mutable internal state such as a cache or interpolator.
    Function local_pdf = pdf;

    // Thread-local basis evaluators and sampling distributions.
    std::vector<LegendreBasis> basis;
    std::vector<std::uniform_real_distribution<double>> dists;
    basis.reserve(dim);
    dists.reserve(dim);
    for (int d = 0; d < dim; ++d) {
      basis.emplace_back(_orders[d]);
      dists.emplace_back(_domain[d].first, _domain[d].second);
    }

    // Run this thread's portion of the Monte Carlo loop.
    auto partial_sums =
        accumulateTrials(local_pdf, basis, dists, gen, trials_this_thread);

    // Normalize before the critical section — only an addition loop
    // remains inside, keeping contention on the shared array minimal.
    const double inv_total = 1.0 / static_cast<double>(_number_of_trials);
    for (double &s : partial_sums)
      s *= inv_total;

#pragma omp critical
    for (int flat = 0; flat < _number_of_co_efficients; ++flat)
      _co_efficients[flat] += partial_sums[flat];
  }
}

} // namespace nontefonte