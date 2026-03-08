/**
 * @file Domain.h
 * @brief Unified domain representation for 1D, 2D, and 3D
 */

#ifndef LOCAL_FET_DOMAIN_H
#define LOCAL_FET_DOMAIN_H

#include <array>
#include <vector>
#include <utility>
#include <random>

namespace local_fet {

    class Domain {
    public:
        Domain(double x_min, double x_max);
        Domain(double x_min, double x_max,
               double y_min, double y_max);
        Domain(double x_min, double x_max,
               double y_min, double y_max,
               double z_min, double z_max);

        int dim() const;

        std::pair<double, double> bounds(int d) const;

        double measure() const;

        template<typename RNG>
        std::vector<double> sample(RNG& rng) const {
            std::vector<double> point(_dim);
            for (int d = 0; d < _dim; ++d) {
                std::uniform_real_distribution<double> dist(
                        _bounds[2*d], _bounds[2*d + 1]);
                point[d] = dist(rng);
            }
            return point;
        }

        bool contains(const std::vector<double>& point) const;

    private:
        int _dim;
        std::array<double, 6> _bounds;
    };

}

#endif