#include "FourierFET.h"
#include <cmath>
#include <stdexcept>

namespace local_fet {

    FourierFET::FourierFET(const Domain& domain, const std::vector<int>& modes)
            : FETBase(domain, modes)
    {
        size_t dim = _orders.size();

        _omega.resize(dim);
        _strides.resize(dim);

        for (size_t d = 0; d < dim; ++d) {
            auto [min,max] = domain.bounds(d);
            double L = max - min;
            _omega[d] = 2.0 * M_PI / L;
        }

        _num_coeffs = 1;
        for (size_t d = 0; d < dim; ++d)
            _num_coeffs *= (1 + 2 * _orders[d]);

        _strides[dim-1] = 1;
        for (int d = (int)dim - 2; d >= 0; --d)
            _strides[d] = _strides[d+1] * (1 + 2 * _orders[d+1]);

        _coefficients.resize(_num_coeffs,0.0);
        _norm_factors.resize(_num_coeffs);

        computeNormalizationFactors();
    }

    void FourierFET::score(const std::vector<double>& point, double weight)
    {
        if (point.size() != _orders.size())
            throw std::invalid_argument("Point dimension mismatch");

        std::vector<std::vector<double>> basis(_orders.size());

        for (size_t d=0; d<_orders.size(); ++d) {
            auto [min,max] = _domain.bounds(d);
            double x = point[d] - min;

            int n = 1 + 2*_orders[d];
            basis[d].resize(n);

            basis[d][0] = 1.0;

            for (int k=1; k<=_orders[d]; ++k) {
                basis[d][2*k-1] = cos(k*_omega[d]*x);
                basis[d][2*k]   = sin(k*_omega[d]*x);
            }
        }

        for (size_t idx=0; idx<_num_coeffs; ++idx) {

            auto mi = fourierMultiIndex(idx);

            double val = 1.0;
            for (size_t d=0; d<_orders.size(); ++d)
                val *= basis[d][mi[d]];

            _coefficients[idx] += weight * val;
        }
    }

    double FourierFET::reconstruct(const std::vector<double>& point) const
    {
        std::vector<std::vector<double>> basis(_orders.size());

        for (size_t d=0; d<_orders.size(); ++d) {
            auto [min,max] = _domain.bounds(d);
            double x = point[d] - min;

            int n = 1 + 2*_orders[d];
            basis[d].resize(n);

            basis[d][0] = 1.0;

            for (int k=1; k<=_orders[d]; ++k) {
                basis[d][2*k-1] = cos(k*_omega[d]*x);
                basis[d][2*k]   = sin(k*_omega[d]*x);
            }
        }

        double result = 0.0;

        for (size_t idx=0; idx<_num_coeffs; ++idx) {

            auto mi = fourierMultiIndex(idx);

            double val = 1.0;
            for (size_t d=0; d<_orders.size(); ++d)
                val *= basis[d][mi[d]];

            result += _coefficients[idx] * val;
        }

        return result;
    }

    void FourierFET::computeNormalizationFactors()
    {
        double measure = _domain.measure();

        for (size_t idx=0; idx<_num_coeffs; ++idx) {

            auto mi = fourierMultiIndex(idx);

            double norm = 1.0;

            for (size_t d=0; d<_orders.size(); ++d)
                norm *= (mi[d]==0) ? 1.0 : 2.0;

            _norm_factors[idx] = norm / measure;
        }
    }

    std::vector<int> FourierFET::fourierMultiIndex(size_t flat) const
    {
        std::vector<int> mi(_orders.size());

        for (size_t d=0; d<_orders.size(); ++d)
            mi[d] = (flat / _strides[d]) % (1 + 2*_orders[d]);

        return mi;
    }

}