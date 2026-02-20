#include "Domain.h"

namespace nontefonte {

Point3D cross_product(const Point3D &a, const Point3D &b) {
  return {a._y * b._z - a._z * b._y, a._z * b._x - a._x * b._z,
          a._x * b._y - a._y * b._x};
}

double dot_product(const Point3D &a, const Point3D &b) { return (a * b).sum(); }

double get_distance(const Point3D &a, const Point3D &b) {
  auto d = a - b;
  return std::sqrt((d * d).sum());
}

std::ostream &operator<<(std::ostream &os, const Point3D &p) {
  os << p._x << ", " << p._y << ", " << p._z;
  return os;
}

Domain::Domain() : _is_problem_2D(false) {}

Domain::Domain(std::vector<Point3D> points, bool is_2D)
    : _is_problem_2D(is_2D), _domain(std::move(points)) {
  if (_domain.size() < 3)
    throw std::invalid_argument("Domain: at least 3 vertices are required.");
  setBoundingBox();
}

Domain::Domain(std::vector<Point2D> points) : _is_problem_2D(true) {
  if (points.size() < 3)
    throw std::invalid_argument("Domain: at least 3 vertices are required.");
  _domain.assign(points.begin(), points.end()); // upcast Point2D -> Point3D
  setBoundingBox();
}

void Domain::setDomain(std::vector<Point3D> points) {
  if (points.size() < 3)
    throw std::invalid_argument("Domain: at least 3 vertices are required.");
  _domain = std::move(points);
  setBoundingBox();
}

void Domain::setDomain(std::vector<Point2D> points) {
  if (!_is_problem_2D)
    throw std::logic_error("Domain: cannot set 2D points on a 3D domain.");
  if (points.size() < 3)
    throw std::invalid_argument("Domain: at least 3 vertices are required.");
  _domain.assign(points.begin(), points.end());
  setBoundingBox();
}

void Domain::setFaces(std::vector<Face3D> faces) {
  if (_is_problem_2D)
    throw std::logic_error("Domain: faces are not used in 2D mode.");
  if (faces.empty())
    throw std::invalid_argument("Domain: face list must not be empty.");
  _faces = std::move(faces);
}

void Domain::setBoundingBox() {
  if (_domain.empty())
    throw std::logic_error(
        "Domain: cannot compute bounding box — no vertices.");

  constexpr double INF = std::numeric_limits<double>::max();
  _bbox = {INF, -INF, INF, -INF, INF, -INF};

  for (const auto &p : _domain) {
    _bbox.xMin = std::min(_bbox.xMin, p.x());
    _bbox.xMax = std::max(_bbox.xMax, p.x());
    _bbox.yMin = std::min(_bbox.yMin, p.y());
    _bbox.yMax = std::max(_bbox.yMax, p.y());
    if (!_is_problem_2D) {
      _bbox.zMin = std::min(_bbox.zMin, p.z());
      _bbox.zMax = std::max(_bbox.zMax, p.z());
    }
  }
  if (_is_problem_2D)
    _bbox.zMin = _bbox.zMax = 0.0;
}

bool Domain::isInsidePolygon2D(const Point3D &p,
                               const std::vector<Point3D> &poly) const {
  int winding = 0;
  const std::size_t n = poly.size();

  for (std::size_t i = 0; i < n; ++i) {
    const Point3D &a = poly[i];
    const Point3D &b = poly[(i + 1) % n];

    if (a.y() <= p.y()) {
      if (b.y() > p.y()) { // upward crossing
        double cross = (b.x() - a.x()) * (p.y() - a.y()) -
                       (p.x() - a.x()) * (b.y() - a.y());
        if (cross > 0.0)
          ++winding;
      }
    } else {
      if (b.y() <= p.y()) { // downward crossing
        double cross = (b.x() - a.x()) * (p.y() - a.y()) -
                       (p.x() - a.x()) * (b.y() - a.y());
        if (cross < 0.0)
          --winding;
      }
    }
  }
  return winding != 0;
}

std::optional<double> Domain::rayTriangleIntersect(const Point3D &orig,
                                                   const Point3D &dir,
                                                   const Point3D &v0,
                                                   const Point3D &v1,
                                                   const Point3D &v2) const {
  constexpr double EPS = 1e-9;

  Point3D e1 = v1 - v0;
  Point3D e2 = v2 - v0;
  Point3D h = cross_product(dir, e2);
  double det = dot_product(e1, h);

  if (std::abs(det) < EPS)
    return std::nullopt; // ray parallel to triangle

  double invDet = 1.0 / det;
  Point3D s = orig - v0;
  double u = dot_product(s, h) * invDet;
  if (u < 0.0 || u > 1.0)
    return std::nullopt;

  Point3D q = cross_product(s, e1);
  double v = dot_product(dir, q) * invDet;
  if (v < 0.0 || u + v > 1.0)
    return std::nullopt;

  double t = dot_product(e2, q) * invDet;
  if (t < EPS)
    return std::nullopt; // intersection behind origin

  return t;
}

bool Domain::isInsideMesh3D(const Point3D &p, const std::vector<Point3D> &verts,
                            const std::vector<Face3D> &faces) const {
  const Point3D rayDir{1.0, 0.0, 0.0};
  int hits = 0;
  for (const auto &f : faces)
    if (rayTriangleIntersect(p, rayDir, verts[f.i0], verts[f.i1], verts[f.i2]))
      ++hits;
  return (hits % 2) == 1;
}

Point3D Domain::samplePointDomain(std::size_t maxAttempts) {
  if (!_is_problem_2D && _faces.empty())
    throw std::logic_error(
        "samplePointDomain: 3D domain has no faces — call setFaces() first.");

  std::uniform_real_distribution<double> distX(_bbox.xMin, _bbox.xMax);
  std::uniform_real_distribution<double> distY(_bbox.yMin, _bbox.yMax);
  std::uniform_real_distribution<double> distZ(_bbox.zMin, _bbox.zMax);

  for (std::size_t attempt = 0; attempt < maxAttempts; ++attempt) {
    Point3D candidate{distX(_rng), distY(_rng),
                      _is_problem_2D ? 0.0 : distZ(_rng)};

    bool inside = _is_problem_2D ? isInsidePolygon2D(candidate, _domain)
                                 : isInsideMesh3D(candidate, _domain, _faces);
    if (inside)
      return candidate;
  }

  throw std::runtime_error("samplePointDomain: exceeded max attempts — domain "
                           "may be too thin or degenerate. Please increase "
                           "maxAttempts more than 100000");
}

Point3D Domain::scaleSampledPointToLegendreDomain(const Point3D &point) const {
  auto mapAxis = [](double val, double lo, double hi) -> double {
    if (hi <= lo)
      throw std::domain_error("scaleSampledPointToLegendreDomain: degenerate "
                              "bounding box axis (min >= max).");
    return 2.0 * (val - lo) / (hi - lo) - 1.0;
  };

  /* for now this is not what I want. I will need to discuss this with Paul
   * before implementing this. But the over all idea is we would draw a stencil
   * by keeping the sampled point in the center. We would then determine the
   * at which points the stencil intersects the domain. Based on the cord length
   * we would scale the sampled point
   */
  double lx = mapAxis(point.x(), _bbox.xMin, _bbox.xMax);
  double ly = mapAxis(point.y(), _bbox.yMin, _bbox.yMax);
  double lz = _is_problem_2D ? 0.0 : mapAxis(point.z(), _bbox.zMin, _bbox.zMax);

  return {lx, ly, lz};
}

} // namespace nontefonte