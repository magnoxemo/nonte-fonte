#include "FETEngine.h"
#include "PDF.h"


#include <iostream>
#include <omp.h>
#include <random>

nontefonte::FETEngine::FETEngine(std::vector<int> orders, long number_of_trials)
    : _orders(orders), _number_of_trials(number_of_trials) {
  for (size_t i = 0; i < _orders.size(); ++i)
    _domain.emplace_back(-1.0, 1.0);
  _domain.shrink_to_fit();
  initCoefficients();
}

nontefonte::FETEngine::FETEngine(std::vector<int> orders, long number_of_trials,
                                 std::vector<std::pair<double, double>> domain)
    : _orders(orders), _number_of_trials(number_of_trials), _domain(domain) {
  _domain.shrink_to_fit();
  initCoefficients();
}

void nontefonte::FETEngine::initCoefficients() {
  int total = 1;
  for (int o : _orders)
    total *= (o + 1);
  _number_of_co_efficients = total;
  _co_efficients.assign(total, 0.0);
}

std::vector<int> nontefonte::FETEngine::multiIndex(int flat) const {
  int dim = _orders.size();
  std::vector<int> idx(dim);
  for (int d = dim - 1; d >= 0; --d) {
    idx[d] = flat % (_orders[d] + 1);
    flat /= (_orders[d] + 1);
  }
  return idx;
}

std::vector<double> nontefonte::FETEngine::samplePoint(
    std::vector<std::uniform_real_distribution<double>> &dists,
    std::default_random_engine &gen) const {
  std::vector<double> x(_orders.size());
  for (size_t d = 0; d < x.size(); ++d)
    x[d] = dists[d](gen);
  return x;
}

double nontefonte::FETEngine::evaluatePDF(Function &pdf,
                                          const std::vector<double> &x) const {
  switch (x.size()) {
  case 1:
    return pdf(x[0]);
  case 2:
    return pdf(x[0], x[1]);
  case 3:
    return pdf(x[0], x[1], x[2]);
  default:
    throw std::runtime_error("Dimension >3 not supported");
  }
}

double
nontefonte::FETEngine::basisProduct(const std::vector<LegendreBasis> &basis,
                                    const std::vector<double> &x,
                                    const std::vector<int> &mi) const {
  double score = 1.0;
  for (size_t d = 0; d < basis.size(); ++d)
    score *= basis[d](static_cast<double >(x[d]), mi[d]);
  return score;
}

std::vector<double> nontefonte::FETEngine::accumulateTrials(
    Function &pdf, std::vector<LegendreBasis> &basis,
    std::vector<std::uniform_real_distribution<double>> &dists,
    std::default_random_engine &gen) const {
  std::vector<double> sums(_number_of_co_efficients, 0.0);

#pragma omp for
  for (long i = 0; i < _number_of_trials; ++i) {
    auto x = samplePoint(dists, gen);
    double pdf_val = evaluatePDF(pdf, x);

    for (int flat = 0; flat < _number_of_co_efficients; ++flat) {
      auto mi = multiIndex(flat);
      double score = basisProduct(basis, x, mi);
      sums[flat] += score * pdf_val / _number_of_trials;
    }
  }
  return sums;
}

void nontefonte::FETEngine::runSimulation(Function &pdf) {
  int dim = _orders.size();

#pragma omp parallel
  {
    std::default_random_engine gen(std::random_device{}() +
                                   omp_get_thread_num());
    auto local_pdf = pdf;

    std::vector<LegendreBasis> basis;
    std::vector<std::uniform_real_distribution<double>> dists;
    for (int d = 0; d < dim; ++d) {
      basis.emplace_back(_orders[d]);
      dists.emplace_back(_domain[d].first, _domain[d].second);
    }

    auto thread_sums = accumulateTrials(local_pdf, basis, dists, gen);

#pragma omp critical
    for (int flat = 0; flat < _number_of_co_efficients; ++flat)
      _co_efficients[flat] += thread_sums[flat];
  }
}