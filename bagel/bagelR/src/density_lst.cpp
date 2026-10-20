
#include "density_lst.h"


real_type density_lst(const real_type& x, const real_type& nu, const real_type& mu,const real_type& sigma)
{
  boost::math::students_t dist(nu);
  return boost::math::pdf(dist,(x-mu)/sigma)/sigma;
}



