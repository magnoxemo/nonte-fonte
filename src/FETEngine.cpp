#include "FETEngine.h"
#include "PDF.h"

#include <random>
#include <omp.h>

nontefonte::FETEngine::FETEngine(std::vector<int> orders, long number_of_trials):_orders(orders),_number_of_trials(number_of_trials){

  int total_order = 1;
  for (size_t index = 0; index < _orders.size(); ++index) {
    total_order =total_order* (_orders[index] + 1);
    if (_domain.empty())
      _domain.emplace_back(-1.0, 1.0);
  }
  _number_of_co_efficients = total_order;
  _co_efficients.assign(total_order, 0.0);
  _domain.shrink_to_fit();
}

nontefonte::FETEngine::FETEngine(std::vector<int> orders, long number_of_trials,
                                 std::vector<std::pair<double, double>> domain)
    : _orders(orders), _number_of_trials(number_of_trials), _domain(domain) {

  int total_order = 1;
  for (size_t index = 0; index < _orders.size(); ++index) {
    total_order =total_order* (_orders[index] + 1);
  }
  _number_of_co_efficients = total_order;
  _co_efficients.resize(_number_of_co_efficients, 0);
  _domain.shrink_to_fit();
}

void nontefonte::FETEngine::runSimulation(Function &pdf) {

#pragma omp declare reduction(                                                 \
        vector_plus : std::vector<double> : std::transform(                    \
                omp_out.begin(), omp_out.end(), omp_in.begin(),                \
                    omp_out.begin(), std::plus<double>()))                     \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


  int dim = pdf.getVariableSize();
  if (dim != _orders.size())
    throw std::runtime_error("Dimension mismatch");

  std::vector<LegendreBasis> basis;
  std::vector<std::uniform_real_distribution<double>> distributions;

  for (int d = 0; d < dim; ++d) {
    basis.emplace_back(LegendreBasis(_orders[d]));
    distributions.emplace_back(_domain[d].first, _domain[d].second);
  }

  std::vector<double> global_coeff(_number_of_co_efficients, 0.0);

#pragma omp parallel
  {
    std::vector<double> local_coeff(_number_of_co_efficients, 0.0);
    std::default_random_engine generator(std::random_device{}() +
                                         omp_get_thread_num());

    for (long i = 0; i < _number_of_trials; ++i) {
      std::vector<double> x(dim);

      for (int d = 0; d < dim; ++d)
        x[d] = distributions[d](generator);

      double pdf_value;

      if (dim == 1)
        pdf_value = pdf(x[0]);
      else if (dim == 2)
        pdf_value = pdf(x[0], x[1]);
      else if (dim == 3)
        pdf_value = pdf(x[0], x[1], x[2]);
      else
        throw std::runtime_error("Dimension >3 not supported");

#pragma omp parallel for reduction (vector_plus:_co_efficients)
      for (int _i = 0; _i < _orders[0]; _i++) {
        for (int _j = 0; _j < _orders[1]; _j++) {
          for (int _k = 0; _k < _orders[2]; _k++) {
            double score = 0;
            if (dim == 1)
              score = basis[0](x[0], _i);
            else if (dim == 2)
              score = basis[0](x[0], _i) *
                      basis[1](x[1], _j);
            else if (dim == 3)
              score = basis[0](x[0], _i) *
                      basis[1](x[1], _j) *
                      basis[2](x[2], _k);
            else
              throw std::runtime_error("Dimension >3 not supported");

            _co_efficients[_i * _orders[0] + _j * _orders[1] + _k] += score * pdf_value ;
          }
        }
      }
    }
  }
}
