#include "Domain.h"
#include <stdexcept>

namespace local_fet {

    Domain::Domain(double x_min, double x_max)
            : _dim(1), _bounds{{x_min, x_max, 0, 0, 0, 0}}
    {
      if (x_min >= x_max) {
        throw std::invalid_argument("Invalid 1D domain");
      }
    }

    Domain::Domain(double x_min, double x_max,
                   double y_min, double y_max)
            : _dim(2), _bounds{{x_min, x_max, y_min, y_max, 0, 0}}
    {
      if (x_min >= x_max || y_min >= y_max) {
        throw std::invalid_argument("Invalid 2D domain");
      }
    }

    Domain::Domain(double x_min, double x_max,
                   double y_min, double y_max,
                   double z_min, double z_max)
            : _dim(3), _bounds{{x_min, x_max, y_min, y_max, z_min, z_max}}
    {
      if (x_min >= x_max || y_min >= y_max || z_min >= z_max) {
        throw std::invalid_argument("Invalid 3D domain");
      }
    }

    int Domain::dim() const {
      return _dim;
    }

    std::pair<double, double> Domain::bounds(int d) const {
      return {_bounds[2*d], _bounds[2*d + 1]};
    }

    double Domain::measure() const {
      double m = 1.0;
      for (int d = 0; d < _dim; ++d) {
        m *= (_bounds[2*d + 1] - _bounds[2*d]);
      }
      return m;
    }

    bool Domain::contains(const std::vector<double>& point) const {
      if (static_cast<int>(point.size()) != _dim)
        return false;

      for (int d = 0; d < _dim; ++d) {
        if (point[d] < _bounds[2*d] ||
            point[d] > _bounds[2*d + 1]) {
          return false;
        }
      }
      return true;
    }

}