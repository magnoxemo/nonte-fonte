/**
 * @file TallyBase.h
 * @brief Abstract base class for all tally types
 */

#ifndef LOCAL_FET_TALLY_BASE_H
#define LOCAL_FET_TALLY_BASE_H

#include "Domain.h"
#include <vector>

namespace nonte_fonte {

/**
 * @class TallyBase
 * @brief Abstract interface for all tallies
 *
 * Each tally knows:
 * - Its own domain
 * - How to check if a point is inside
 * - How to score a weighted sample
 * - How to finalize after all samples
 */
    class TallyBase {
    public:
        /**
         * @brief Construct a tally with a domain
         */
        explicit TallyBase(const Domain& domain) : _domain(domain) {}

        /// Virtual destructor
        virtual ~TallyBase() = default;

        /**
         * @brief Check if point is in this tally's domain
         */
        virtual bool contains(const std::vector<double>& point) const {
            return _domain.contains(point);
        }

        /**
         * @brief Score a weighted sample
         *
         * @param point Sample point (dimensionality matches domain)
         * @param weight Sample weight (typically pdf(point) * domain_measure)
         */
        virtual void score(const std::vector<double>& point, double weight) = 0;

        /**
         * @brief Finalize tally after all samples
         *
         * @param total_weight Sum of all weights across all samples
         */
        virtual void finalize(double total_weight) = 0;

        /// Get tally's domain
        const Domain& domain() const { return _domain; }

    protected:
        Domain _domain;  ///< Tally's spatial domain
    };

} // namespace local_fet

#endif // LOCAL_FET_TALLY_BASE_H