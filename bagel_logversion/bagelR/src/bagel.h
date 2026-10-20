#ifndef ___BAGEL_H___
#define ___BAGEL_H___

#include "probability_type.h"
#include "time_type.h"
#include "real_type.h"
#include "ratio_type.h"
#include "tau_type.h"
#include "time_type.h"
#include "matrix_type.h"
#include "particle_type.h"
#include "KL_divergence_type.h"
#include "particle_type.h"
#include "theorem_2.h"
#include "theorem_3.h"
#include "theorem_1.h"
#include "theorem_4.h"
#include "normal_density.h"
#include <list>
#include <cmath>
#include <limits>
#include <algorithm>
#include "prune.h"

#include "model_type.h"



template<typename noise, KL_divergence_type KL_divergence>
struct bagel_type
{

  time_type t;
  probability_type p;
  probability_type p0;
  real_type s;
  int max_num_particles; 
  
  particle_type<noise,KL_divergence> initial_particle;
  std::list<particle_type<noise,KL_divergence> > particles;


  bagel_type(const model_type& model,
	     const noise& noise_structure,
	     const probability_type& p0,
	     const probability_type& p,
	     const int& n)
  {
    t = 1;
    max_num_particles = n;
    initial_particle = particle_type<noise,KL_divergence>(model,noise_structure,0,p0,p);
  }


  // no-change posterior weight on the log scale; retained internally so it
  // stays finite and informative even when exp() of it would underflow to zero
  real_type log_weight_0_t() const
  {
    return particles.front().log_weight;
  }

  real_type weight_0_t() const
  {
    return std::exp(particles.front().log_weight);
  }

  std::list<std::list<ratio_type> > ratios() const
  {
    std::list<std::list<ratio_type> >  particle_ratios;
    for(auto& p : particles)
    {
      particle_ratios.push_back(p.ratios);
    }
    return particle_ratios;
  }


  bagel_type& update(const real_type& y)
  {
    
    if(t == 1)
      {
	initial_particle = theorem_2(initial_particle,t);
	initial_particle = theorem_1(initial_particle,initial_particle,t);
	initial_particle = theorem_4(initial_particle,t,y);
	initial_particle = theorem_3(initial_particle,t,y);
	initial_particle.log_weight = std::log(initial_particle.p0);
	particles.push_back(initial_particle);
      }
    else
      {
	// new particle
	initial_particle = particles.front();
	particles.push_back(theorem_2(initial_particle,t));
	for(auto& p : particles)
	  {
	    p = theorem_1(initial_particle,p,t);
	  }

	for(auto& p : particles)
	  {
	    p = theorem_4(p,t,y);
	    p = theorem_3(p,t,y);
	  }
	// normalise on the log scale (log-sum-exp) so no intermediate
	// exponentiation can underflow before the weights are normalised
	real_type max_log_weight = -std::numeric_limits<real_type>::infinity();
	for(auto& p : particles)
	  {
	    max_log_weight = std::max(max_log_weight,p.log_weight);
	  }
	real_type sum_of_scaled_weights = 0.0;
	for(auto& p : particles)
	  {
	    sum_of_scaled_weights += std::exp(p.log_weight - max_log_weight);
	  }
	real_type log_sum_of_weights = max_log_weight + std::log(sum_of_scaled_weights);
	for(auto& p : particles)
	  {
	    p.log_weight = p.log_weight - log_sum_of_weights;
	  }

	particles = prune(particles,max_num_particles);
      }
    t = t + 1; 
    return *this;
  }
};

template<typename noise, KL_divergence_type KL_divergence>
bagel_type<noise,KL_divergence>& update(bagel_type<noise,KL_divergence>& bagel, const real_type& y)
{
  return bagel.update(y);
}

template<typename noise, KL_divergence_type KL_divergence>
real_type weight_0_t(const bagel_type<noise,KL_divergence>& bagel)
{
  return bagel.weight_0_t();
}

template<typename noise, KL_divergence_type KL_divergence>
real_type log_weight_0_t(const bagel_type<noise,KL_divergence>& bagel)
{
  return bagel.log_weight_0_t();
}

template<typename noise, KL_divergence_type KL_divergence>
std::list<std::list<ratio_type> > ratios(const bagel_type<noise,KL_divergence>& bagel)
{
  return bagel.ratios();
}


#endif
