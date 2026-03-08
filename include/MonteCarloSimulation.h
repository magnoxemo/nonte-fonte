/**
 * @file MCSimulation.h
 * @brief Unified Monte Carlo simulation for N-dimensional PDFs
 */

#ifndef LOCAL_FET_MC_SIMULATION_H
#define LOCAL_FET_MC_SIMULATION_H

#include "TallyBase.h"
#include "Domain.h"

#include <vector>
#include <functional>
#include <iostream>

namespace nonte_fonte {

    void printLogo()
    {
        std::cout <<R"(
        ._  __      ._.      .-----.          .-.
        | \ | | ___ | |_ ___ |  ___|__  _ __  | |_ ___
        |  \| |/ _ \| __/ _ \| |_ / _ \| '_ \ | __/ _ \
        | |\  | (_) | ||  __/|  _| (_) | | | || ||  __/
        |_| \_|\___/ \__\___||_|  \___/|_| |_| \__\___|

                     n o t e - f o n t e
        A FET testing mini app.

        Ebny Walid Ahammed
        Computational Nuclear Engineering Research Lab
        Dept of Nuclear Engineering and Engineering Physics
        University of Wisconsin-Madison.
        )" << std::endl;
    }

    class MCSimulation {
    public:
        MCSimulation(const nonte_fonte::Domain& global_domain,
                     std::function<double(const std::vector<double>&)> target_pdf,
                     long n_samples,
                     unsigned seed = 0);

        void addTally(nonte_fonte::TallyBase* tally);

        void run();

        long numSamples() const { return _n_samples; }

        int dim() const { return _domain.dim(); }

    private:
        nonte_fonte::Domain _domain;
        std::function<double(const std::vector<double>&)> _pdf;
        long _n_samples;
        unsigned _master_seed;
        std::vector<TallyBase*> _tallies;
    };

}

#endif