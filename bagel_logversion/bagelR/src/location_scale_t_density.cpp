
#include "location_scale_t_density.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>

real_type location_scale_t_density(const real_type& x, const real_type& nu, const real_type& mu,const real_type& sigma)
{
  boost::math::students_t dist(nu);
  return boost::math::pdf(dist,(x-mu)/sigma)/sigma;
}

// log of location_scale_t_density, evaluated directly on the log scale so the
// predictive weight update never underflows for extreme observations
real_type log_location_scale_t_density(const real_type& x, const real_type& nu, const real_type& mu,const real_type& sigma)
{
  real_type z = (x - mu)/sigma;
  return std::lgamma(0.5*(nu + 1.0)) - std::lgamma(0.5*nu)
       - 0.5*std::log(nu*boost::math::constants::pi<real_type>())
       - 0.5*(nu + 1.0)*std::log1p(z*z/nu)
       - std::log(sigma);
}
