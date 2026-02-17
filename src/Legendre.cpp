#include "Legendre.h"

nontefonte::LegendreBasis::LegendreBasis(const int n, Dimension dimension ):_N(n),_dim(dimension){
    if (_N<0)
        throw std::runtime_error("Number of order must be non negative");
    if (_dim != Dimension::x || _dim != Dimension::y ||_dim != Dimension::z ){
        throw std::runtime_error("invalid dimension");
    }
};


template<typename T>
double nontefonte::LegendreBasis::operator()(T x, const int order)  {
    if (order == 0) {
        return 1.0;
    }
    if (order == 1) {
        return x;
    }

    double p_prev2 = 1.0;
    double p_prev1 = x;
    double p_curr = 0.0;

    for (int i = 1; i < order; ++i) {

        p_curr = ((2.0 * i + 1.0) * x * p_prev1 - i * p_prev2) / (i + 1.0);
        p_prev2 = p_prev1;
        p_prev1 = p_curr;
    }

    return p_curr;
}