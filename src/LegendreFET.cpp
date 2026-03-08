#include "LegendreFET.h"
#include <cmath>

namespace nonte_fonte {

    LegendreFET::LegendreFET(const Domain& domain, const std::vector<int>& orders)
            : FETBase(domain, orders) {
        computeNormalizationFactors();
    }

    void LegendreFET::score(const std::vector<double>& point, double weight)
    {
        std::vector<double> xi = mapToBasisDomain(point);

        std::vector<std::vector<double>> P(_orders.size());

        for (size_t d = 0; d < _orders.size(); ++d) {
            P[d].resize(_orders[d] + 1);
            for (int n = 0; n <= _orders[d]; ++n) {
                P[d][n] = evaluateLegendre(xi[d], n);
            }
        }

        for (int idx = 0; idx < _num_coeffs; ++idx) {

            std::vector<int> mi = multiIndex(idx);

            double basis_val = 1.0;

            for (size_t d = 0; d < _orders.size(); ++d) {
                basis_val *= P[d][mi[d]];
            }

            _coefficients[idx] += weight * basis_val;
        }
    }

    double LegendreFET::reconstruct(const std::vector<double>& point) const
    {
        std::vector<double> xi = mapToBasisDomain(point);

        std::vector<std::vector<double>> P(_orders.size());

        for (size_t d = 0; d < _orders.size(); ++d) {
            P[d].resize(_orders[d] + 1);
            for (int n = 0; n <= _orders[d]; ++n) {
                P[d][n] = evaluateLegendre(xi[d], n);
            }
        }

        double result = 0.0;

        for (int idx = 0; idx < _num_coeffs; ++idx) {

            std::vector<int> mi = multiIndex(idx);

            double basis_val = 1.0;

            for (size_t d = 0; d < _orders.size(); ++d) {
                basis_val *= P[d][mi[d]];
            }

            result += _coefficients[idx] * basis_val;
        }

        return result;
    }

    void LegendreFET::computeNormalizationFactors()
    {
        _norm_factors.resize(_num_coeffs);

        for (int idx = 0; idx < _num_coeffs; ++idx) {

            std::vector<int> mi = multiIndex(idx);

            double norm = 1.0;

            for (size_t d = 0; d < _orders.size(); ++d) {

                auto bounds = _domain.bounds(d);

                double L_d = bounds.second - bounds.first;

                norm *= (2.0 * mi[d] + 1.0) / L_d;
            }

            _norm_factors[idx] = norm;
        }
    }

    double LegendreFET::evaluateLegendre(double x, int n) const
    {
        if (n == 0) return 1.0;
        if (n == 1) return x;

        double p_prev2 = 1.0;
        double p_prev1 = x;
        double p_curr = 0.0;

        for (int k = 1; k < n; ++k) {
            p_curr = ((2.0 * k + 1.0) * x * p_prev1 - k * p_prev2) / (k + 1.0);
            p_prev2 = p_prev1;
            p_prev1 = p_curr;
        }

        return p_curr;
    }

}