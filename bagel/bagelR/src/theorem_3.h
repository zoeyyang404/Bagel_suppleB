

#ifndef ___THEOREM_3_H___
#define ___THEOREM_3_H___

#include "particle_type.h"
#include "time_type.h"
#include "real_type.h"


template <typename noise_type, KL_divergence_type KL_divergence>
requires requires { requires std::same_as<noise_type,known_variance>; }
particle_type<noise_type,KL_divergence>& theorem_3_noise(particle_type<noise_type,KL_divergence>& particle_t, const matrix& e, const matrix& Q)
{
  return particle_t;
}


template <typename noise_type, KL_divergence_type KL_divergence>
requires requires { requires std::same_as<noise_type,unknown_variance>; }
particle_type<noise_type,KL_divergence>& theorem_3_noise(particle_type<noise_type,KL_divergence>& particle_t, const matrix& e, const matrix& Q)
{
  particle_t.noise_structure.nu = particle_t.noise_structure.nu + 0.5;
  particle_t.noise_structure.iota = particle_t.noise_structure.iota + 0.5*e(0,0)*e(0,0)/Q(0,0);
  return particle_t;

}


// modifies exiting particles - mutating
template<typename noise, KL_divergence_type KL_divergence>
particle_type<noise,KL_divergence>& theorem_3(particle_type<noise,KL_divergence>& particle_t, const time_type& t, const real_type& val)
{
  matrix y(1,1);
  y(0,0) = val;
  prior_type prior_t = particle_t.model.prior_function(t);
  matrix mu = particle_t.post.mu;
  matrix sigma = particle_t.post.sigma;
  tau_type tau = particle_t.tau;
  matrix H_t = particle_t.model.feature_vector_function(t,tau);
  matrix I = matrix::Identity(1,1);
  matrix e = y - H_t.transpose() * mu;  
  matrix Q = H_t.transpose() * sigma * H_t + I;
  real_type rQ = Q(0,0);
  matrix A = sigma * H_t / rQ;
  sigma = sigma - A * A.transpose() * rQ;
  mu = mu + A * e;
  particle_t.post.mu = mu;
  particle_t.post.sigma = sigma;
  particle_t = theorem_3_noise(particle_t,e,Q);
  return particle_t;
}





#endif

