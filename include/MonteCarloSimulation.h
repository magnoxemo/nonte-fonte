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
         _      ____  _      _____ _____     _____ ____  _      _____ _____
        / \  /|/  _ \/ \  /|/__ __Y  __/    / ___//  _ \/ \  /|/__ __Y  __/
        | |\ ||| / \|| |\ ||  / \ |  \  --- |  __\| / \|| |\ ||  / \ |  \
        | | \||| \_/|| | \||  | | |  /_ --- | |   | \_/|| | \||  | | |  /_
        \_/  \|\____/\_/  \|  \_/ \____\    \_/   \____/\_/  \|  \_/ \____\

                     n o t e - f o n t e
        A FET testing mini app.

        Ebny Walid Ahammed
        Computational Nuclear Engineering Research Group
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

        std::vector<TallyBase*> getTallies(){
            return _tallies;
        }
    private:
        nonte_fonte::Domain _domain;
        std::function<double(const std::vector<double>&)> _pdf;
        long _n_samples;
        unsigned _master_seed;
        std::vector<TallyBase*> _tallies;
    };

}

#endif
