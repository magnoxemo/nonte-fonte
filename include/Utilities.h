#ifndef NONTE_FONTE_UTILITIES_H
#define NONTE_FONTE_UTILITIES_H

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>

namespace nonte_fonte {

    struct Statistics {
        double mean = 0.0;
        double variance = 0.0;
        double std_dev = 0.0;
        int n_samples = 0;

        void update(double value) {
            n_samples++;
            double delta = value - mean;
            mean += delta / n_samples;
            double delta2 = value - mean;
            variance += delta * delta2;
        }

        void finalize() {
            if (n_samples > 1) {
                variance /= (n_samples - 1);
                std_dev = std::sqrt(variance);
            }
        }

        double relative_error() const {
            return (mean != 0.0) ? std_dev / std::abs(mean) : 0.0;
        }
    };

    class CoefficientStatistics {
    private:
        std::vector<Statistics> _coeff_stats;

    public:
        explicit CoefficientStatistics(size_t n_coeffs)
                : _coeff_stats(n_coeffs) {}

        void addSample(const std::vector<double>& coeffs) {
            if (coeffs.size() != _coeff_stats.size()) {
                throw std::runtime_error("Coefficient size mismatch");
            }
            for (size_t i = 0; i < coeffs.size(); ++i) {
                _coeff_stats[i].update(coeffs[i]);
            }
        }

        void finalize() {
            for (auto& stat : _coeff_stats) {
                stat.finalize();
            }
        }

        const std::vector<Statistics>& getStats() const {
            return _coeff_stats;
        }

        std::vector<double> getMeans() const {
            std::vector<double> means(_coeff_stats.size());
            for (size_t i = 0; i < _coeff_stats.size(); ++i) {
                means[i] = _coeff_stats[i].mean;
            }
            return means;
        }

        std::vector<double> getStdDevs() const {
            std::vector<double> stds(_coeff_stats.size());
            for (size_t i = 0; i < _coeff_stats.size(); ++i) {
                stds[i] = _coeff_stats[i].std_dev;
            }
            return stds;
        }

        void printSummary(const std::string& prefix = "") const {
            std::cout << prefix << "Coefficient Statistics (mean ± std_dev):\n";
            for (size_t i = 0; i < _coeff_stats.size(); ++i) {
                const auto& stat = _coeff_stats[i];
                std::cout << prefix << "  c[" << std::setw(3) << i << "]: "
                          << std::setw(12) << std::setprecision(6) << std::fixed << stat.mean
                          << " ± " << std::setw(10) << stat.std_dev;

                if (stat.mean != 0.0) {
                    double rel_err = stat.relative_error() * 100.0;
                    std::cout << "  (" << std::setw(6) << std::setprecision(2) << rel_err << "%)";
                }
                std::cout << "\n";
            }
        }
    };


    inline void printCoefficientsWithUncertainty(
            const CoefficientStatistics& stats,
            const std::string& label = "Coefficients")
    {
        std::cout << "\n" << label << ":\n";
        std::cout << std::setw(8) << "Index"
                  << std::setw(15) << "Mean"
                  << std::setw(15) << "Std Dev"
                  << std::setw(12) << "Rel Err (%)\n";
        std::cout << std::string(50, '-') << "\n";

        const auto& coeff_stats = stats.getStats();
        for (size_t i = 0; i < coeff_stats.size(); ++i) {
            const auto& stat = coeff_stats[i];
            std::cout << std::setw(8) << i
                      << std::setw(15) << std::setprecision(8) << std::fixed << stat.mean
                      << std::setw(15) << std::setprecision(8) << stat.std_dev;

            if (stat.mean != 0.0) {
                double rel_err = stat.relative_error() * 100.0;
                std::cout << std::setw(12) << std::setprecision(4) << rel_err;
            } else {
                std::cout << std::setw(12) << "N/A";
            }
            std::cout << "\n";
        }
    }

} // namespace nonte_fonte

#endif // NONTE_FONTE_UTILITIES_