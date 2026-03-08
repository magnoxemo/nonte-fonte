/**
 * @file ZernikeFET.h
 * @brief Zernike polynomial FET for circular domains
 */

#ifndef LOCAL_FET_ZERNIKE_FET_H
#define LOCAL_FET_ZERNIKE_FET_H

#include "TallyBase.h"
#include <vector>
#include <utility>

namespace nonte_fonte {

    class ZernikeFET : public TallyBase {
    public:
        ZernikeFET(const Domain& domain, double radius, int max_n);

        void score(const std::vector<double>& point, double weight) override;

        void finalize(double total_weight) override;

        double reconstruct(const std::vector<double>& point) const;

        const std::vector<double>& coefficients() const { return _coefficients; }

        const std::vector<std::pair<int,int>>& modes() const { return _modes; }

    private:
        double evaluateZernike(double rho, double theta, int n, int m) const;

        double evaluateRadial(double rho, int n, int m) const;

        double factorial(int n) const;

        void computeNormalization();

        double _radius;
        int _max_n;

        std::vector<std::pair<int,int>> _modes;
        std::vector<double> _coefficients;
        std::vector<double> _norm_factors;
    };

}

#endif