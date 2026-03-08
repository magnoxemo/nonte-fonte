#include "FETBase.h"
#include <stdexcept>

namespace local_fet {

    FETBase::FETBase(const Domain& domain, const std::vector<int>& orders)
            : TallyBase(domain), _orders(orders)
    {
      if (static_cast<int>(orders.size()) != domain.dim()) {
        throw std::invalid_argument(
                "Orders size must match domain dimensionality");
      }

      for (int order : orders) {
        if (order < 0) {
          throw std::invalid_argument(
                  "Orders must be non-negative");
        }
      }

      _num_coeffs = 1;

      for (int order : orders) {
        _num_coeffs *= (order + 1);
      }

      _coefficients.resize(_num_coeffs, 0.0);

    }

    void FETBase::finalize(double total_weight)
    {
      if (total_weight <= 0.0) {
        throw std::runtime_error(
                "Total weight must be positive");
      }

      for (int i = 0; i < _num_coeffs; ++i) {
        _coefficients[i] /= total_weight;
        _coefficients[i] *= _norm_factors[i];
      }
    }

    std::vector<double> FETBase::mapToBasisDomain(
            const std::vector<double>& point) const
    {
      std::vector<double> xi(point.size());

      for (size_t d = 0; d < point.size(); ++d) {

        auto bounds = _domain.bounds(d);

        double min = bounds.first;
        double max = bounds.second;

        xi[d] = 2.0 * (point[d] - min) / (max - min) - 1.0;
      }

      return xi;
    }

    std::vector<int> FETBase::multiIndex(int flat) const
    {
      std::vector<int> mi(_orders.size());

      for (int d = static_cast<int>(_orders.size()) - 1; d >= 0; --d) {
        mi[d] = flat % (_orders[d] + 1);
        flat /= (_orders[d] + 1);
      }

      return mi;
    }

    int FETBase::flatIndex(const std::vector<int>& mi) const
    {
      int idx = 0;
      int stride = 1;

      for (int d = static_cast<int>(_orders.size()) - 1; d >= 0; --d) {
        idx += mi[d] * stride;
        stride *= (_orders[d] + 1);
      }

      return idx;
    }

}