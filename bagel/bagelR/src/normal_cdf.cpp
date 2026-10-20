#include "normal_cdf.h"

#include <iostream>
#include <cmath> 


real_type normal_cdf(const real_type& x)
{
    return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}
