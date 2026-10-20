
#ifndef ___THEOREM_2_H___
#define ___THEOREM_2_H___

#include "particle_type.h"
#include "time_type.h"


// creates new particles - non mutating
template<typename noise, KL_divergence_type KL_divergence>
particle_type<noise,KL_divergence> theorem_2(const particle_type<noise,KL_divergence>& particle_0, const time_type& t)
{

  particle_type<noise,KL_divergence> particle_t = particle_0;

  if(t == 1)
    {
      particle_t.post = particle_0.model.prior_function(t);
      particle_t.tau= t - 1;
      return particle_t;
    }

  // Make everyhting very explicit for now - can take some shortcuts (i.e. no temporary objects) once it is tested.
  
  prior_type prior_1 = particle_0.model.prior_function(1);
  prior_type prior_t = particle_0.model.prior_function(t);
  
  int d1 = prior_1.mu.rows();
  int d2 = prior_t.mu.rows() - d1;
  int n = d1 + d2;

  // resize particle post matrices
  particle_t.post.mu = Eigen::MatrixXd::Zero(n,1);
  particle_t.post.sigma = Eigen::MatrixXd::Zero(n,n);
  

  
  matrix sigma_beta_beta_t = prior_t.sigma.block(0,0,d1,d1);     
  matrix sigma_beta_gamma_t = prior_t.sigma.block(0,d1,d1,d2);


  
  matrix sigma_gamma_beta_t = prior_t.sigma.block(d1,0,d2,d1);
  matrix sigma_gamma_gamma_t = prior_t.sigma.block(d1,d1,d2,d2);
  
  matrix mu_beta_t = prior_t.mu.block(0,0,d1,1);     
  matrix mu_gamma_t = prior_t.mu.block(d1,0,d2,1);

  matrix B = sigma_beta_gamma_t.transpose() * sigma_beta_beta_t.inverse();
  matrix top_left = particle_0.post.sigma.block(0,0,d1,d1);
  matrix top_right = B * particle_0.post.sigma.block(0,0,d1,d1);
  matrix bottom_left = top_right.transpose();
  matrix bottom_right = sigma_gamma_gamma_t + B * (particle_0.post.sigma.block(0,0,d1,d1) - sigma_beta_beta_t) * B.transpose();
  matrix top = particle_0.post.mu.block(0,0,d1,1);
  matrix bottom = mu_gamma_t + B * (particle_0.post.mu.block(0,0,d1,1) - mu_beta_t);

  particle_t.post.mu = matrix::Zero(d1 + d2,1);
  particle_t.post.sigma = matrix::Zero(d1 + d2,d1 + d2);
  
  particle_t.post.sigma.block(0,0,d1,d1) = top_left;    
  particle_t.post.sigma.block(0,d1,d1,d2) = top_right.transpose();
  particle_t.post.sigma.block(d1,0,d2,d1) = bottom_left.transpose();
  particle_t.post.sigma.block(d1,d1,d2,d2) = bottom_right;

  particle_t.post.mu.block(0,0,d1,1) = top;
  particle_t.post.mu.block(d1,0,d2,1) = bottom;

  particle_t.tau= t - 1;
  particle_t.noise_structure = particle_0.noise_structure;


  return particle_t;
  
}



#endif
