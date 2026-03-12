#ifndef LOCAL_FET_TALLY_BASE_H
#define LOCAL_FET_TALLY_BASE_H

#include "Domain.h"
#include <vector>

namespace nonte_fonte {

    class TallyBase {
    public:

        explicit TallyBase(const Domain& domain) : _domain(domain) {}

        virtual ~TallyBase() = default;

        virtual bool contains(const std::vector<double>& point) const {
            return _domain.contains(point);
        }
        virtual void score(const std::vector<double>& point, double weight) = 0;
        virtual void finalize(double total_weight) = 0;
        const Domain& domain() const { return _domain; }

    protected:
        Domain _domain;
    };

} // namespace local_fet

#endif // LOCAL_FET_TALLY_BASE_H