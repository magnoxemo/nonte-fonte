#include <omp.h>
#include <random>
#include <vector>
#include <iostream>

#include "PDF.h"
#include "Legendre.h"



#pragma omp declare reduction(vector_plus : std::vector<double> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
        initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


int main(int argc, char* argv[]){

    const long int N = std::stoi(argv[1]);
    const int order = std::stoi(argv[2]);
    float average = 0.0f;
    nontefonte::Function pdf("x*x-3*x",{"x"});
    nontefonte::LegendreBasis legendre_polynomials(order, nontefonte::Dimension::x);

    std::vector<double> co_efficients (order+1, 0.0);
    #pragma omp parallel
    {
        // creating a local PDF and random number generator
        // otherwise race condition complains
        auto local_pdf = pdf;
        auto local_legendre_polynomials = legendre_polynomials;

        std::default_random_engine random_number_generator( std::random_device{}() + omp_get_thread_num());
        std::uniform_real_distribution<double> distribution(-1.0f, 1.0f);

        auto rng = [&random_number_generator, & distribution](){
            return distribution(random_number_generator);
        };

        #pragma omp parallel for reduction( vector_plus: co_efficients)
        for (long i = 0; i < N; i++) {
            auto x = rng();
            auto pdf_value = local_pdf(x);
            for (int j =0 ; j< order; j++){
                co_efficients[j] += pdf_value* local_legendre_polynomials(x, j);
            }
        }
    }

    printf("[");
    for (int i =0; i<order; i++) {
        std::cout<<  (co_efficients[i]/ N)<<", ";
    }
    printf("]\n\n");

    return 0;
}
