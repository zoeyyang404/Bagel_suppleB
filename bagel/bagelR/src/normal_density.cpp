
#include "normal_density.h"
#include <boost/math/distributions/normal.hpp>

real_type normal_density(const real_type& x,const real_type& mu,const real_type& sd)
{
  boost::math::normal Z; // standard normal distribution
  return boost::math::pdf(Z,(x - mu)/sd)/sd;
}


