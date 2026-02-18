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
    public:
        LegendreBasis(const int n);


        template <typename T>
        double operator() (T x, const int order) {
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
        };

        const int _N;
    };
}

#endif //NONTE_FONTE_LEGENDRE_H
