#include "Legendre.h"

nontefonte::LegendreBasis::LegendreBasis(const int n, Dimension dimension ):_N(n),_dim(dimension){
    if (_N<0)
        throw std::runtime_error("Number of order must be non negative");
//    if (_dim != Dimension::x || _dim != Dimension::y || _dim != Dimension::z ){
//        throw std::runtime_error("invalid dimension");
//    }
};


