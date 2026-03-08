/**
 * @file LegendreFET.h
 * @brief Legendre polynomial FET for 1D, 2D, and 3D
 */

#ifndef LOCAL_FET_LEGENDRE_FET_H
#define LOCAL_FET_LEGENDRE_FET_H

#include "FETBase.h"
#include <vector>

namespace local_fet {

    class LegendreFET : public FETBase {
    public:
        LegendreFET(const Domain& domain, const std::vector<int>& orders);

        void score(const std::vector<double>& point, double weight) override;

        double reconstruct(const std::vector<double>& point) const override;

    protected:
        void computeNormalizationFactors() override;

    private:
        double evaluateLegendre(double x, int n) const;
    };

}

#endif