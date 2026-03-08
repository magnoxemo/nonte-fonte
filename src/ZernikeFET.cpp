#include "ZernikeFET.h"
#include <cmath>
#include <stdexcept>

namespace nonte_fonte {

    ZernikeFET::ZernikeFET(const Domain& domain, double radius, int max_n)
            : TallyBase(domain), _radius(radius), _max_n(max_n)
    {
        if (domain.dim() != 2)
            throw std::invalid_argument("Domain must be 2D");

        if (radius <= 0.0 || max_n < 0)
            throw std::invalid_argument("Invalid radius or order");

        for (int n = 0; n <= max_n; ++n) {
            for (int m = -n; m <= n; ++m) {
                if ((n - std::abs(m)) % 2 == 0)
                    _modes.push_back({n, m});
            }
        }

        _coefficients.resize(_modes.size(), 0.0);

        computeNormalization();
    }

    void ZernikeFET::score(
            const std::vector<double>& point,
            double weight)
    {
        double x = point[0];
        double y = point[1];

        double r = std::sqrt(x*x + y*y);

        if (r > _radius)
            return;

        double theta = std::atan2(y, x);
        double rho = r / _radius;

        for (size_t i = 0; i < _modes.size(); ++i) {

            auto mode = _modes[i];

            double Z = evaluateZernike(
                    rho,
                    theta,
                    mode.first,
                    mode.second
            );

            _coefficients[i] += weight * Z;
        }
    }

    void ZernikeFET::finalize(double total_weight)
    {
        if (total_weight <= 0.0)
            throw std::runtime_error("Total weight must be positive");

        for (size_t i = 0; i < _coefficients.size(); ++i) {
            _coefficients[i] /= total_weight;
            _coefficients[i] *= _norm_factors[i];
        }
    }

    double ZernikeFET::reconstruct(
            const std::vector<double>& point) const
    {
        double x = point[0];
        double y = point[1];

        double r = std::sqrt(x*x + y*y);

        if (r > _radius)
            return 0.0;

        double theta = std::atan2(y, x);
        double rho = r / _radius;

        double result = 0.0;

        for (size_t i = 0; i < _modes.size(); ++i) {

            auto mode = _modes[i];

            double Z = evaluateZernike(
                    rho,
                    theta,
                    mode.first,
                    mode.second
            );

            result += _coefficients[i] * Z;
        }

        return result;
    }

    double ZernikeFET::evaluateZernike(
            double rho,
            double theta,
            int n,
            int m) const
    {
        double R = evaluateRadial(rho, n, std::abs(m));

        if (m >= 0)
            return R * std::cos(m * theta);
        else
            return R * std::sin(std::abs(m) * theta);
    }

    double ZernikeFET::evaluateRadial(
            double rho,
            int n,
            int m) const
    {
        if ((n - m) % 2 != 0)
            return 0.0;

        double result = 0.0;

        int k_max = (n - m) / 2;

        for (int k = 0; k <= k_max; ++k) {

            double sign = (k % 2 == 0) ? 1.0 : -1.0;

            double coeff =
                    sign * factorial(n - k) /
                    ( factorial(k)
                      * factorial((n + m)/2 - k)
                      * factorial((n - m)/2 - k) );

            result += coeff * std::pow(rho, n - 2*k);
        }

        return result;
    }

    double ZernikeFET::factorial(int n) const
    {
        static std::vector<double> cache = {1.0};

        while (static_cast<int>(cache.size()) <= n)
            cache.push_back(cache.back() * cache.size());

        return (n >= 0 && n < static_cast<int>(cache.size()))
               ? cache[n]
               : 1.0;
    }

    void ZernikeFET::computeNormalization()
    {
        _norm_factors.resize(_modes.size());

        double area = M_PI * _radius * _radius;

        for (size_t i = 0; i < _modes.size(); ++i) {
            int n = _modes[i].first;
            _norm_factors[i] = (n + 1.0) / area;
        }
    }

}