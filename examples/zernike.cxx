#include "Domain.h"
#include "TallyBase.h"
#include "MonteCarloSimulation.h"
#include "FETBase.h"
#include "ZernikeFET.h"

#include <iostream>
#include <fstream>
#include <cmath>

double gaussian_circular(const std::vector<double>& x) {
    double r_sq = x[0]*x[0] + x[1]*x[1];
    return std::exp(-r_sq) / M_PI;
}


void write_output(std::ofstream& file,double radius, const nonte_fonte::ZernikeFET& basis){

    file << "x,y,true_pdf,legendre_pdf\n";

    for (double x = -radius; x <= radius; x += 0.01) {
        for (double y = -radius; y <= radius; y += 0.01) {
            std::vector<double> point = {x, y};
            double r = std::sqrt(x*x + y*y);
            if (r <= radius) {
                file << x << "," << y << ","
                     << gaussian_circular(point) << ","
                     << basis.reconstruct(point) << "\n";
            }
        }
    }
    file.close();
}

int main() {

    double radius = 1;
    long n_samples = 5000000;
    nonte_fonte::Domain domain(-radius, radius, -radius, radius);
    nonte_fonte::ZernikeFET basis(domain, radius, 10);

    std::cout << "Number of Zernike modes: " << basis.modes().size() << "\n\n";

    nonte_fonte::MCSimulation mc(domain, gaussian_circular, n_samples);
    mc.addTally(&basis);
    mc.run();

    std::cout << "\n";

    std::ofstream file("zernike_results.csv");

    write_output(file, radius, basis);


    return 0;
}