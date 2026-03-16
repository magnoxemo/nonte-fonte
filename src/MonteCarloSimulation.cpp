#include "MonteCarloSimulation.h"
#include "TallyBase.h"
#include "Domain.h"

#include <random>
#include <iostream>
#include <stdexcept>
#include <omp.h>

namespace nonte_fonte {

    MCSimulation::MCSimulation(const nonte_fonte::Domain& global_domain,
                               std::function<double(const std::vector<double>&)> target_pdf,
                               long n_samples,
                               unsigned seed)
            : _domain(global_domain),
              _pdf(target_pdf),
              _n_samples(n_samples),
              _master_seed(seed == 0 ? std::random_device{}() : seed)
    {
      if (_n_samples <= 0) {
        throw std::invalid_argument("Number of samples must be positive");
      }
    }

    void MCSimulation::addTally(nonte_fonte::TallyBase* tally)
    {
      if (tally == nullptr) {
        throw std::invalid_argument("Cannot add null tally");
      }

      _tallies.push_back(tally);
    }

    void MCSimulation::run()
    {
      if (_tallies.empty()) {
        std::cout << "Warning: No tallies registered\n";
        return;
      }

      const int dim = _domain.dim();
      const double measure = _domain.measure();

      nonte_fonte::printLogo();
      std::cout<<"\n\n";
      std::cout << "=================================================\n";
      std::cout << "             Monte Carlo Simulation\n";
      std::cout << "=================================================\n";
      std::cout << "     Dimensionality: " << dim << "D\n";
      std::cout << "     Domain measure: " << measure << "\n";
      std::cout << "     Total samples: " << _n_samples << "\n";
      std::cout << "     Tallies: " << _tallies.size() << "\n";

      int num_threads = omp_get_max_threads();

      std::cout << "    OpenMP threads: " << num_threads << "\n";
      std::cout << "========================================\n\n";

      std::vector<double> thread_weights(num_threads, 0.0);

      #pragma omp parallel
      {
        int thread_id = omp_get_thread_num();

        std::mt19937 rng(_master_seed + thread_id);

        double local_weight = 0.0;

        long samples_per_thread = _n_samples / num_threads;
        long start = thread_id * samples_per_thread;
        long end = (thread_id == num_threads - 1) ? _n_samples : start + samples_per_thread;

        for (long i = start; i < end; ++i) {

          std::vector<double> point = _domain.sample(rng);
          double weight = _pdf(point) * measure;

          local_weight += weight;

          for (auto* tally : _tallies) {
            if (tally->contains(point)) {
              #pragma omp critical
              {
                tally->score(point, weight);
              }
            }
          }

        }

        thread_weights[thread_id] = local_weight;
      }

      double total_weight = 0.0;

      for (double w : thread_weights) {
        total_weight += w;
      }

      std::cout << "\nTotal weight: " << total_weight << "\n";
      std::cout << "Finalizing tallies...\n";

      for (auto* tally : _tallies) {
        tally->finalize(total_weight);
      }

      std::cout << "Simulation complete.\n";
      std::cout << "========================================\n";
    }

}