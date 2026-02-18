#include "FETEngine.h"
#include "PDF.h"

#include <random>
#include <omp.h>
#include <iostream>

nontefonte::FETEngine::FETEngine(std::vector<int> orders, long number_of_trials):_orders(orders),_number_of_trials(number_of_trials){

  int total_order = 1;

  for (size_t index = 0; index < _orders.size(); ++index) {
    total_order = total_order* (_orders[index] + 1);
    _domain.emplace_back(-1.0, 1.0);
  }

  std::cout<<total_order;
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

void nontefonte::FETEngine::runSimulation(Function & pdf) {


#pragma omp parallel
  {
    std::default_random_engine generator(std::random_device{}() +
                                         omp_get_thread_num());

    int dim = _orders.size();
    auto local_pdf = pdf ;

    std::vector<LegendreBasis> basis;
    std::vector<std::uniform_real_distribution<double>> distributions;

    for (int d = 0; d < dim; ++d) {
      basis.emplace_back(LegendreBasis(_orders[d]));
      distributions.emplace_back(_domain[d].first, _domain[d].second);
    }

#pragma omp for
    for (long i = 0; i < _number_of_trials; ++i) {
      std::vector<double> x(local_pdf.getVariableSize(),0);

      for (int d = 0; d < dim; ++d)
        x[d] = distributions[d](generator);

      double pdf_value;

      if (local_pdf.getVariableSize() == 1)
        pdf_value = local_pdf(x[0]);
      else if (local_pdf.getVariableSize() == 2)
        pdf_value = local_pdf(x[0], x[1]);
      else if (local_pdf.getVariableSize() == 3)
        pdf_value = local_pdf(x[0], x[1], x[2]);
      else
        throw std::runtime_error("local_pdf.getVariableSize() >3 not supported");


#pragma omp critical
      {
      for (int _i = 0; _i < _orders[0]; _i++) {
        for (int _j = 0; _j < _orders[1]; _j++) {
          for (int _k = 0; _k < _orders[2]; _k++) {
            double score = 0;
            int idx ;
            if (dim == 1) {
              score = basis[0](x[0], _i);
              idx = _i;
            }
            else if (dim == 2) {
              score = basis[0](x[0], _i) * basis[1](x[1], _j);
              idx = _i * (_orders[1] + 1) + _j;
            }
            else if (dim == 3) {
              score = basis[0](x[0], _i) * basis[1](x[1], _j) * basis[2](x[2], _k);
              int idx = _i * ( _orders[1] + 1 * _orders[2] + 1) + _j * _orders[2] + 1 + _k;

            }
            else
              throw std::runtime_error("Dimension >3 not supported");

            _co_efficients[idx] += score * pdf_value / _number_of_trials;
          }
        }
      }
    }
    }
  }
}
