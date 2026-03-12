#ifndef LOCAL_FET_HISTOGRAM_TALLY_H
#define LOCAL_FET_HISTOGRAM_TALLY_H

#include "TallyBase.h"
#include <vector>
#include <array>
#include <utility>

namespace nonte_fonte {

    class HistogramTally1D : public TallyBase {
    public:
        HistogramTally1D(const Domain& domain, const std::vector<double>& mesh);

        void score(const std::vector<double>& point, double weight) override;

        void finalize(double total_weight) override;

        std::vector<double> binCenters() const;

        const std::vector<double>& values() const { return _counts; }

        const std::vector<double>& mesh() const { return _mesh; }

    private:
        std::vector<double> _mesh;
        std::vector<double> _counts;
        std::vector<double> _bin_widths;
    };

    class HistogramTally2D : public TallyBase {
    public:
        HistogramTally2D(const Domain& domain,
                         const std::vector<double>& x_edges,
                         const std::vector<double>& y_edges);

        void score(const std::vector<double>& point, double weight) override;

        void finalize(double total_weight) override;

        double value(size_t ix, size_t iy) const;

        const std::vector<double>& values() const { return _counts; }

        std::pair<size_t, size_t> numBins() const { return {_nx, _ny}; }

    private:
        std::vector<double> _x_edges;
        std::vector<double> _y_edges;
        size_t _nx, _ny;
        std::vector<double> _counts;
        std::vector<double> _areas;
    };

    class HistogramTallyUnstructured : public TallyBase {
    public:
        HistogramTallyUnstructured(
                const Domain& domain,
                const std::vector<std::array<double,2>>& mesh_points);

        void score(const std::vector<double>& point, double weight) override;

        void finalize(double total_weight) override;

        const std::vector<double>& values() const { return _counts; }

        const std::vector<std::array<double,2>>& meshPoints() const { return _mesh_points; }

    private:
        std::vector<std::array<double,2>> _mesh_points;
        std::vector<double> _counts;
        std::vector<double> _areas;
    };

}

#endif