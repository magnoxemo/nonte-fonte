#include "Domain.h"


namespace nontefonte{

    inline Point3D cross_product(const Point3D &point_a, const Point3D &point_b) {
        return {point_a._y * point_b._z - point_a._z * point_b._y,
                point_a._z * point_b._x - point_a._x * point_b._z,
                point_a._x * point_b._y - point_a._y * point_b._x};
    }

    inline double dot_product(const Point3D &point_a, const Point3D &point_b) {
        return (point_a * point_b).sum();
    }

    inline double get_distance(const Point3D &point_a, const Point3D &point_b) {

        auto l = point_a - point_b;
        return sqrt((l * l).sum());
    }

    inline std::ostream &operator<<(std::ostream &os, const Point3D &p) {
        os << p._x << ", " << p._y << ", " << p._z;
        return os;
    }
}