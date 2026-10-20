
#ifndef ___PARTICLE_TYPE_H___
#define ___PARTICLE_TYPE_H___


#include "prior_type.h"
#include "post_type.h"
#include "real_type.h"
#include "tau_type.h"
#include "ratio_type.h"
#include "noise_type.h"
#include "KL_divergence_type.h"
#include <list>

#include "model_type.h"



template<typename noise, KL_divergence_type KL_divergence>
struct particle_type
{
  model_type model;
  noise noise_structure;
  prior_function_type prior_function;
  feature_vector_function_type feature_vector_function;
  tau_type tau;
  real_type p0;
  real_type p;
  post_type post;
  real_type weight;
  std::list<ratio_type> ratios;

  particle_type() {};


  particle_type& operator=(const particle_type& other)
    {
      model = other.model;
      prior_function = other.prior_function;
      feature_vector_function = other.feature_vector_function;
      noise_structure = other.noise_structure;
      tau = other.tau;
      p0 = other.p0;
      p = other.p;
      post = other.post;
      weight = other.weight;
      ratios = other.ratios;
      return *this;
    }

  particle_type(const model_type& _model,
		const noise& _noise_structure,
		const tau_type& _tau,
		const real_type& _p0,
		const real_type& _p) 
  {
    model = _model;
    noise_structure = _noise_structure;
    tau = _tau;
    p0 = _p0;
    p = _p;
    weight = 1.0;
    ratios.push_front(1.0);
  }

  
};












#endif

