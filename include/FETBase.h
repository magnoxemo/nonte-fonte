#ifndef LOCAL_FET_FET_BASE_H
#define LOCAL_FET_FET_BASE_H

#include "TallyBase.h"
#include <vector>

namespace nonte_fonte {

    class FETBase : public TallyBase {
    public:
        FETBase(const Domain& domain, const std::vector<int>& orders);

        const std::vector<int>& orders() const { return _orders; }

        const std::vector<double>& coefficients() const { return _coefficients; }

        std::vector<double>& coefficients() { return _coefficients; }

        void finalize(double total_weight) override;

        virtual double reconstruct(const std::vector<double>& point) const = 0;

    protected:
        virtual void computeNormalizationFactors() = 0;

        std::vector<double> mapToBasisDomain(const std::vector<double>& point) const;

        std::vector<int> multiIndex(int flat) const;

        int flatIndex(const std::vector<int>& mi) const;

        std::vector<int> _orders;
        int _num_coeffs;
        std::vector<double> _coefficients;
        std::vector<double> _norm_factors;
    };

}

#endif