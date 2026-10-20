

#ifndef ___THEOREM_1_H___
#define ___THEOREM_1_H___

#include "particle_type.h"
#include "time_type.h"
#include "real_type.h"


template<typename noise, KL_divergence_type KL_divergence>
particle_type<noise,KL_divergence>& theorem_1(particle_type<noise,KL_divergence>& particle_1,particle_type<noise,KL_divergence>& particle_t, const time_type& t)
{
  real_type weight;
  if(particle_t.tau == 0)
    {
      if(t == 1)
	{
	  particle_t.weight = 1.0;
	}

      if(t == 2)
	{
	  particle_t.weight = particle_t.p0;
	}
      return(particle_t); 
    }
  if(particle_t.tau == 1)
    {
      if(t == 2)
	{
	  particle_t.weight = 1.0 - particle_t.p0;
	  return particle_t;
	}    
    }
    if(particle_t.tau == t - 1)
      {
	particle_t.weight = particle_1.weight*(1.0 - particle_t.p0)/(particle_t.p0 * real_type(t -1));
      }
    else
      {
	particle_t.weight = particle_t.weight*real_type(t-2)/real_type(t-1);
      }
    return particle_t;
}




#endif
