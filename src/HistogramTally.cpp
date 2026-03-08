#include "HistogramTally.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace nonte_fonte {


    HistogramTally1D::HistogramTally1D(
            const Domain& domain,
            const std::vector<double>& mesh)
            : TallyBase(domain), _mesh(mesh){
        if (domain.dim() != 1)
            throw std::invalid_argument("Domain must be 1D");

        if (mesh.size() < 2)
            throw std::invalid_argument("Mesh must have at least 2 points");

        if (!std::is_sorted(mesh.begin(), mesh.end()))
            throw std::invalid_argument("Mesh must be sorted");

        _counts.resize(mesh.size(), 0.0);
        _bin_widths.resize(mesh.size() - 1);

        for (size_t i = 0; i < mesh.size() - 1; ++i)
            _bin_widths[i] = mesh[i+1] - mesh[i];
    }

    void HistogramTally1D::score(
            const std::vector<double>& point,
            double weight)
    {
        double x = point[0];

        auto it = std::lower_bound(_mesh.begin(), _mesh.end(), x);

        if (it == _mesh.end())
            _counts.back() += weight;

        else if (it == _mesh.begin() && x < _mesh[0])
            _counts[0] += weight;

        else {
            size_t idx = std::distance(_mesh.begin(), it);
            if (idx > 0) idx--;
            _counts[idx] += weight;
        }
    }

    void HistogramTally1D::finalize(double total_weight)
    {
        if (total_weight <= 0.0)
            throw std::runtime_error("Total weight must be positive");

        for (size_t i = 0; i < _counts.size(); ++i) {
            _counts[i] /= total_weight;
            if (i < _bin_widths.size())
                _counts[i] /= _bin_widths[i];
        }
    }

    std::vector<double> HistogramTally1D::binCenters() const
    {
        std::vector<double> centers(_mesh.size() - 1);

        for (size_t i = 0; i < centers.size(); ++i)
            centers[i] = 0.5 * (_mesh[i] + _mesh[i+1]);

        return centers;
    }

/* ---------------- 2D Histogram ---------------- */

    HistogramTally2D::HistogramTally2D(
            const Domain& domain,
            const std::vector<double>& x_edges,
            const std::vector<double>& y_edges)
            : TallyBase(domain),
              _x_edges(x_edges),
              _y_edges(y_edges)
    {
        if (domain.dim() != 2)
            throw std::invalid_argument("Domain must be 2D");

        if (x_edges.size() < 2 || y_edges.size() < 2)
            throw std::invalid_argument("Edges must have at least 2 points");

        if (!std::is_sorted(x_edges.begin(), x_edges.end()) ||
            !std::is_sorted(y_edges.begin(), y_edges.end()))
            throw std::invalid_argument("Edges must be sorted");

        _nx = x_edges.size() - 1;
        _ny = y_edges.size() - 1;

        _counts.resize(_nx * _ny, 0.0);
        _areas.resize(_nx * _ny);

        for (size_t iy = 0; iy < _ny; ++iy)
            for (size_t ix = 0; ix < _nx; ++ix) {
                double dx = x_edges[ix+1] - x_edges[ix];
                double dy = y_edges[iy+1] - y_edges[iy];
                _areas[iy*_nx + ix] = dx * dy;
            }
    }

    void HistogramTally2D::score(
            const std::vector<double>& point,
            double weight)
    {
        double x = point[0];
        double y = point[1];

        auto ix_it = std::lower_bound(_x_edges.begin(), _x_edges.end(), x);
        auto iy_it = std::lower_bound(_y_edges.begin(), _y_edges.end(), y);

        size_t ix = (ix_it == _x_edges.begin())
                    ? 0
                    : std::distance(_x_edges.begin(), ix_it) - 1;

        size_t iy = (iy_it == _y_edges.begin())
                    ? 0
                    : std::distance(_y_edges.begin(), iy_it) - 1;

        ix = std::min(ix, _nx - 1);
        iy = std::min(iy, _ny - 1);

        _counts[iy*_nx + ix] += weight;
    }

    void HistogramTally2D::finalize(double total_weight)
    {
        if (total_weight <= 0.0)
            throw std::runtime_error("Total weight must be positive");

        for (size_t i = 0; i < _counts.size(); ++i) {
            _counts[i] /= total_weight;
            _counts[i] /= _areas[i];
        }
    }

    double HistogramTally2D::value(size_t ix, size_t iy) const
    {
        return _counts[iy*_nx + ix];
    }

/* -------- Unstructured Histogram -------- */

    HistogramTallyUnstructured::HistogramTallyUnstructured(
            const Domain& domain,
            const std::vector<std::array<double,2>>& mesh_points)
            : TallyBase(domain),
              _mesh_points(mesh_points)
    {
        if (domain.dim() != 2)
            throw std::invalid_argument("Domain must be 2D");

        if (mesh_points.empty())
            throw std::invalid_argument("Mesh must have at least one point");

        _counts.resize(mesh_points.size(), 0.0);

        _areas.resize(mesh_points.size());

        double total_area = domain.measure();
        double area = total_area / mesh_points.size();

        std::fill(_areas.begin(), _areas.end(), area);
    }

    void HistogramTallyUnstructured::score(
            const std::vector<double>& point,
            double weight)
    {
        double x = point[0];
        double y = point[1];

        size_t nearest = 0;
        double min_dist = std::numeric_limits<double>::max();

        for (size_t i = 0; i < _mesh_points.size(); ++i) {

            double dx = x - _mesh_points[i][0];
            double dy = y - _mesh_points[i][1];
            double dist = dx*dx + dy*dy;

            if (dist < min_dist) {
                min_dist = dist;
                nearest = i;
            }
        }

        _counts[nearest] += weight;
    }

    void HistogramTallyUnstructured::finalize(double total_weight)
    {
        if (total_weight <= 0.0)
            throw std::runtime_error("Total weight must be positive");

        for (size_t i = 0; i < _counts.size(); ++i) {
            _counts[i] /= total_weight;
            _counts[i] /= _areas[i];
        }
    }

}