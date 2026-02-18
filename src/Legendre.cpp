#include "Legendre.h"

nontefonte::LegendreBasis::LegendreBasis(const int n ):_N(n){
    if (_N<0)
        throw std::runtime_error("Number of order must be non negative");
};


