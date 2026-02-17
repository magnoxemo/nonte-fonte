#include "PDF.h"
#include "omp.h"
#include "random"


int main(int argc, char* argv[]){

    const long int N = std::stoi(argv[1]);
    float average = 0.0f;
    nontefonte::Function pdf("sin(x)/cos(y)+(y-x/cos(x))^2 - sqrt(x*x)", {"x", "y"});

    #pragma omp parallel
    {
        // creating a local PDF and random number generator
        // otherwise race condition complains
        auto local_pdf = pdf;
        std::default_random_engine random_number_generator( std::random_device{}() + omp_get_thread_num());
        std::uniform_real_distribution<float> distribution(-1.0f, 1.0f);
        auto rng = [&random_number_generator, & distribution](){
            return distribution(random_number_generator);
        };

        #pragma omp for reduction(+:average)
        for (long i = 0; i < N; i++) {
            average += local_pdf(rng(), rng());
        }
    }


    printf("average = %f\n",average/N);

    return 0;

}
