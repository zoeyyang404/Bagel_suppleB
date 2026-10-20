
#include "density_norm.h"

real_type density_norm(const real_type& x,const real_type& mu,const real_type& sd)
{
  boost::math::normal Z; // standard normal distribution
  return boost::math::pdf(Z,(x - mu)/sd)/sd;
}


