#ifndef LOCAL_FET_FOURIER_FET_H
#define LOCAL_FET_FOURIER_FET_H

#include "FETBase.h"
#include <vector>

namespace local_fet {

    class FourierFET : public FETBase {
    public:
        FourierFET(const Domain& domain, const std::vector<int>& modes);

        void score(const std::vector<double>& point, double weight) override;
        double reconstruct(const std::vector<double>& point) const override;

    protected:
        void computeNormalizationFactors() override;

    private:
        std::vector<int> fourierMultiIndex(size_t flat) const;

        std::vector<double> _omega;
        std::vector<size_t> _strides;
        size_t _num_coeffs;
    };

}

#endif