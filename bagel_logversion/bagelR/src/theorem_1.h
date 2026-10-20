

#ifndef ___THEOREM_1_H___
#define ___THEOREM_1_H___

#include "particle_type.h"
#include "time_type.h"
#include "real_type.h"

#include <cmath>


// weight recursion carried out on the log scale so the no-change weight can
// underflow-free reach arbitrarily small values
template<typename noise, KL_divergence_type KL_divergence>
particle_type<noise,KL_divergence>& theorem_1(particle_type<noise,KL_divergence>& particle_1,particle_type<noise,KL_divergence>& particle_t, const time_type& t)
{
  if(particle_t.tau == 0)
    {
      if(t == 1)
	{
	  particle_t.log_weight = 0.0;
	}

      if(t == 2)
	{
	  particle_t.log_weight = std::log(particle_t.p0);
	}
      return(particle_t);
    }
  if(particle_t.tau == 1)
    {
      if(t == 2)
	{
	  particle_t.log_weight = std::log1p(-particle_t.p0);
	  return particle_t;
	}
    }
    if(particle_t.tau == t - 1)
      {
	particle_t.log_weight = particle_1.log_weight + std::log1p(-particle_t.p0) - std::log(particle_t.p0) - std::log(real_type(t-1));
      }
    else
      {
	particle_t.log_weight = particle_t.log_weight + std::log(real_type(t-2)) - std::log(real_type(t-1));
      }
    return particle_t;
}




#endif
