
#include "location_scale_t_density.h"
#include <boost/math/distributions/students_t.hpp>

real_type location_scale_t_density(const real_type& x, const real_type& nu, const real_type& mu,const real_type& sigma)
{
  boost::math::students_t dist(nu);
  return boost::math::pdf(dist,(x-mu)/sigma)/sigma;
}
