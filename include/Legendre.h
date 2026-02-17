#ifndef NONTE_FONTE_LEGENDRE_H
#define NONTE_FONTE_LEGENDRE_H

#include <stdexcept>


namespace nontefonte{

    enum Dimension{
        x,
        y,
        z
    };

    /* this class should hold a sets of legendre polynomial
     * functions up to order n */
    class LegendreBasis{
        LegendreBasis(const int n, Dimension dimension );
        ~LegendreBasis()=default;

        template <typename T>
        double operator() (T x, const int order);

        const int _N;
        Dimension _dim;
    };
}

#endif //NONTE_FONTE_LEGENDRE_H
