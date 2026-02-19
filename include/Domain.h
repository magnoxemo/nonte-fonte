#ifndef NONTE_FONTE_DOMAIN_H
#define NONTE_FONTE_DOMAIN_H

#include <cmath>
#include <ostream>
#include <vector>

namespace nontefonte {
struct Point3D {
  double _x, _y, _z;
  Point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}
  constexpr static double _tolerance = 10e-5;

  inline double x() const { return _x; }
  inline double y() const { return _y; }
  inline double z() const { return _z; }

  // Overloading the addition operator (+)
  Point3D operator+(const Point3D &other) const {
    return {_x + other._x, _y + other._y, _z + other._z};
  }

  Point3D operator+(double t) const { return {_x + t, _y + t, _z + t}; }

  // Overloading the addition operator (-)
  Point3D operator-(const Point3D &other) const {
    return {_x - other._x, _y - other._y, _z - other._z};
  }

  Point3D operator-(double t) const { return {_x - t, _y - t, _z - t}; }

  // Overloading the addition operator (*)
  Point3D operator*(const Point3D &other) const {
    return {_x * other._x, _y * other._y, _z * other._z};
  }

  Point3D operator*(double t) const { return {_x * t, _y * t, _z * t}; }

  // Overloading the addition operator
  Point3D operator/(const Point3D &other) const {
    return {_x / other._x, _y / other._y, _z / other._z};
  }

  Point3D operator/(double t) const { return {_x / t, _y / t, _z / t}; }

  bool operator==(const Point3D &other) const {
    Point3D p = {_x - other._x, _y - other._y, _z - other._z};
    if ((p * p).sum() < _tolerance)
      return true;
    return false;
  }

  double sum() const { return _x + _y + _z; }
  double norm() const {
    auto norm = std::sqrt(this->_x * this->_x + this->_y * this->_y +
                          this->_z * this->_z);
    return norm;
  }
};

struct Point2D : Point3D {
        Point2D(double x, double y) : Point3D(x, y, 0) {}
    };


inline Point3D cross_product(const Point3D &point_a, const Point3D &point_b);
inline double dot_product(const Point3D &point_a, const Point3D &point_b);
inline double get_distance(const Point3D &point_a, const Point3D &point_b);
inline std::ostream &operator<<(std::ostream &os, const Point3D &p);


class Domain{
    
    Domain();
};
} // namespace nontefonte

#endif // NONTE_FONTE_DOMAIN_H
