#ifndef NONTE_FONTE_DOMAIN_H
#define NONTE_FONTE_DOMAIN_H

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <ostream>
#include <random>
#include <stdexcept>
#include <variant>
#include <vector>

namespace nontefonte {

// ═══════════════════════════════════════════════════════════════════════════
//  Point types
// ═══════════════════════════════════════════════════════════════════════════

struct Point3D {
  double _x, _y, _z;

  Point3D(double x, double y, double z) : _x(x), _y(y), _z(z) {}

  constexpr static double _tolerance = 1e-9;

  inline double x() const { return _x; }

  inline double y() const { return _y; }

  inline double z() const { return _z; }

  Point3D operator+(const Point3D &other) const {
    return {_x + other._x, _y + other._y, _z + other._z};
  }

  Point3D operator+(double t) const { return {_x + t, _y + t, _z + t}; }

  Point3D operator-(const Point3D &other) const {
    return {_x - other._x, _y - other._y, _z - other._z};
  }

  Point3D operator-(double t) const { return {_x - t, _y - t, _z - t}; }

  Point3D operator*(const Point3D &other) const {
    return {_x * other._x, _y * other._y, _z * other._z};
  }

  Point3D operator*(double t) const { return {_x * t, _y * t, _z * t}; }

  Point3D operator/(const Point3D &other) const {
    return {_x / other._x, _y / other._y, _z / other._z};
  }

  Point3D operator/(double t) const { return {_x / t, _y / t, _z / t}; }

  bool operator==(const Point3D &other) const {
    Point3D p = {_x - other._x, _y - other._y, _z - other._z};
    return (p * p).sum() < _tolerance;
  }

  double sum() const { return _x + _y + _z; }

  double norm() const { return std::sqrt(_x * _x + _y * _y + _z * _z); }
};

struct Point2D : Point3D {
  Point2D(double x, double y) : Point3D(x, y, 0) {}
};

// ═══════════════════════════════════════════════════════════════════════════
//  Supporting types
// ═══════════════════════════════════════════════════════════════════════════

inline Point3D cross_product(const Point3D &a, const Point3D &b);
inline double dot_product(const Point3D &a, const Point3D &b);
inline double get_distance(const Point3D &a, const Point3D &b);
inline std::ostream &operator<<(std::ostream &os, const Point3D &p);


/**
 * @brief A triangular face defined by three indices into a vertex array.
 *
 * Used by the 3D inside test (ray casting). Indices refer to positions in the
 * vector passed to @ref Domain::setDomain(std::vector<Point3D>).
 */
struct Face3D {
  std::size_t i0, ///< Index of the first vertex.
      i1,         ///< Index of the second vertex.
      i2;         ///< Index of the third vertex.
};

/**
 * @brief Axis-aligned bounding box (AABB) for a domain.
 *
 * For 2D domains the z extents are always 0.
 */
struct BoundingBox {
  double xMin = 0.0, xMax = 0.0;
  double yMin = 0.0, yMax = 0.0;
  double zMin = 0.0, zMax = 0.0; ///< Unused for 2D domains.
};

// ═══════════════════════════════════════════════════════════════════════════
//  Domain
// ═══════════════════════════════════════════════════════════════════════════

/**
 * @brief Geometric domain defined by a 2D polygon or a 3D triangulated surface
 * mesh.
 *
 * Internally all points are stored as Point3D (since Point2D inherits from it).
 * The flag @ref _is_problem_2D distinguishes the two cases and gates all
 * dimension-specific behaviour.
 *
 * **2D mode** — vertices define a simple, possibly non-convex polygon
 * (ordered CCW or CW, first vertex not repeated).  Inside test: winding number.
 *
 * **3D mode** — vertices + triangular faces define a *closed* (watertight)
 * surface mesh.  Inside test: Möller–Trumbore ray casting (odd-hit rule).
 *
 * Point sampling uses rejection sampling: candidates are drawn uniformly from
 * the axis-aligned bounding box and accepted only when they pass the inside
 * test. Accepted points can be linearly mapped to the Legendre reference domain
 * [-1,1]^n.
 *
 * @note Dimensionality is fixed at construction and cannot be changed.
 *
 * @example
 * @code
 *   // 2D — non-convex polygon
 *   std::vector<Point3D> verts = {{0,0,0},{2,0,0},{2,1,0},{1,0.5,0},{0,1,0}};
 *   Domain d(verts, true);
 *   Point3D s = d.samplePointDomain();
 *   Point3D l = d.scaleSampledPointToLegendreDomain(s);
 * @endcode
 */
class Domain {
public:
  // ── Constructors ──────────────────────────────────────────────────────────

  /**
   * @brief Default constructor. Creates an empty 3D domain.
   *
   * Call @ref setDomain and @ref setFaces before sampling.
   */
  Domain();

  /**
   * @brief Constructs a domain from a vector of Point3D vertices.
   *
   * Computes the bounding box immediately.
   *
   * @param points      Ordered polygon vertices (2D) or mesh vertices (3D).
   * @param is_2D       Pass @c true for a 2D polygon, @c false for a 3D mesh.
   * @throws std::invalid_argument if fewer than 3 points are provided.
   */
  Domain(std::vector<Point3D> points, bool is_2D);

  /**
   * @brief Convenience constructor for a 2D domain from Point2D vertices.
   *
   * @param points  Ordered 2D polygon vertices (CCW or CW, not closed).
   * @throws std::invalid_argument if fewer than 3 points are provided.
   */
  explicit Domain(std::vector<Point2D> points);

  // ── Setters ───────────────────────────────────────────────────────────────

  /**
   * @brief Replaces the domain vertices and recomputes the bounding box.
   *
   * @param points  New vertex list.
   * @throws std::invalid_argument if fewer than 3 points are provided.
   */
  void setDomain(std::vector<Point3D> points);

  /**
   * @brief Replaces the domain vertices with 2D points and recomputes the
   * bounding box.
   *
   * @param points  New 2D polygon vertices.
   * @throws std::logic_error      if the domain was constructed as 3D.
   * @throws std::invalid_argument if fewer than 3 points are provided.
   */
  void setDomain(std::vector<Point2D> points);

  /**
   * @brief Sets the triangular faces used by the 3D inside test.
   *
   * Must be called before @ref samplePointDomain on a 3D domain.
   * Each face contains three indices (i0, i1, i2) into the vertex array.
   * The mesh must be closed (watertight) for the ray-casting test to be
   * correct.
   *
   * @param faces  Triangle connectivity list.
   * @throws std::logic_error      if the domain is 2D (faces are not used in
   * 2D).
   * @throws std::invalid_argument if @p faces is empty.
   */
  void setFaces(std::vector<Face3D> faces);

  /**
   * @brief Computes and caches the axis-aligned bounding box from stored
   * vertices.
   *
   * Called automatically by constructors and setters. For 2D domains the z
   * extents remain 0.
   *
   * @throws std::logic_error if the vertex list is empty.
   */
  void setBoundingBox();

  // ── Getters ───────────────────────────────────────────────────────────────

  /** @brief Returns the cached bounding box. */
  BoundingBox getBoundingBox() const { return _bbox; }

  /** @brief Returns @c true if the domain is 2D. */
  bool is2D() const { return _is_problem_2D; }

  // ── Sampling ──────────────────────────────────────────────────────────────

  /**
   * @brief Samples a uniformly random point strictly inside the domain.
   *
   * Draws candidates uniformly from the bounding box and rejects those that
   * lie outside the domain, repeating until an interior point is found.
   *
   * - **2D**: inside test is the winding number (handles non-convex polygons).
   * - **3D**: inside test is Möller–Trumbore ray casting with the odd-hit rule
   *           (requires a watertight mesh set via @ref setFaces).
   *
   * @param maxAttempts  Maximum number of candidates before giving up.
   *                     Default is 100 000, which is sufficient for most
   * domains.
   * @return A Point3D (or Point2D-compatible) interior point.
   * @throws std::runtime_error if @p maxAttempts is exceeded.
   * @throws std::logic_error   if a 3D domain has no faces set.
   */
  Point3D samplePointDomain(std::size_t maxAttempts = 100000);

  /**
   * @brief Maps a domain point linearly to the Legendre reference domain [-1,
   * 1]^n.
   *
   * Each axis is transformed independently:
   * @f[
   *   \xi = \frac{2(x - x_{\min})}{x_{\max} - x_{\min}} - 1
   * @f]
   * For 2D domains the z component of the returned point is always 0.
   *
   * @param point  A point expected to lie within the bounding box.
   * @return The mapped point in [-1, 1]^n.
   * @throws std::domain_error if any active bounding box axis is degenerate
   *         (min == max), which would produce a division by zero.
   */
  Point3D scaleSampledPointToLegendreDomain(const Point3D &point) const;

protected:
  /// @brief True if the domain is 2D (polygon), false if 3D (mesh).
  const bool _is_problem_2D;

  /// @brief Stored domain vertices. Point2D instances are up-cast to Point3D.
  std::vector<Point3D> _domain;

  /// @brief Triangle face connectivity for 3D domains. Empty for 2D domains.
  std::vector<Face3D> _faces;

  /// @brief Cached axis-aligned bounding box.
  BoundingBox _bbox;

  /// @brief Mersenne Twister RNG, seeded from hardware entropy at construction.
  std::mt19937 _rng{std::random_device{}()};

  // ── Internal geometry helpers ─────────────────────────────────────────────

  /**
   * @brief Winding-number point-in-polygon test (2D).
   *
   * Counts how many times the polygon boundary winds around @p p. A non-zero
   * winding number means the point is inside. Works correctly for non-convex
   * and multiply-wound polygons; immune to the ray-through-vertex degeneracy
   * that afflicts simple ray-casting methods.
   *
   * Only the x and y components of all points are used.
   *
   * @param p     Candidate point.
   * @param poly  Ordered polygon vertices (the domain's _domain vector).
   * @return @c true if @p p is strictly inside the polygon.
   */
  bool isInsidePolygon2D(const Point3D &p,
                         const std::vector<Point3D> &poly) const;

  /**
   * @brief Ray-casting point-in-mesh test (3D).
   *
   * Fires a ray from @p p in the +X direction and counts how many triangles
   * it intersects. An odd count indicates the point is inside the closed mesh
   * (Jordan curve theorem generalised to 3D). Requires a watertight mesh.
   *
   * @param p      Candidate point.
   * @param verts  Mesh vertices.
   * @param faces  Triangle connectivity.
   * @return @c true if @p p is strictly inside the mesh.
   */
  bool isInsideMesh3D(const Point3D &p, const std::vector<Point3D> &verts,
                      const std::vector<Face3D> &faces) const;

  /**
   * @brief Möller–Trumbore ray–triangle intersection.
   *
   * Tests whether the ray (orig + t·dir, t > 0) intersects the triangle
   * (v0, v1, v2) and returns the distance @c t if so.
   *
   * @param orig      Ray origin.
   * @param dir       Ray direction (need not be unit length).
   * @param v0,v1,v2  Triangle vertices.
   * @return          Distance @c t > 0 to the intersection, or @c std::nullopt.
   */
  std::optional<double> rayTriangleIntersect(const Point3D &orig,
                                             const Point3D &dir,
                                             const Point3D &v0,
                                             const Point3D &v1,
                                             const Point3D &v2) const;
};
}
#endif // namespace nontefonte