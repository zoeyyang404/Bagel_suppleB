
#ifndef ___PRUNE_H___
#define ___PRUNE_H___


#include "particle_type.h"
#include "pruning_type.h"
#include "KL_divergence_type.h"
#include "normal_cdf.h"
#include <list>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <algorithm>
#include <iterator>
template <typename T>
std::tuple<matrix,matrix> transform(const T& a)
{
  matrix mu = a.post.mu;
  matrix sigma = a.post.sigma;
  time_type tau = a.tau;
  matrix A = a.model.transformer_function(tau);
  mu = A*mu;
  sigma = A * sigma * A.transpose();
  return std::make_tuple(mu,sigma);
}

inline double gaussian_kl(const matrix& mu_i, const matrix& Sig_i,
                          const matrix& mu_j, const matrix& Sig_j,
                          const real_type& sigma)
{
  matrix I = matrix::Identity(Sig_i.rows(), Sig_i.cols());
  matrix Sig_j_inv = Sig_j.inverse();
  double result = (Sig_j_inv * Sig_i - I).trace();
  result = result + (1.0/(sigma*sigma)) * ((mu_j - mu_i).transpose() * Sig_j_inv * (mu_j - mu_i))(0,0);
  result = result + std::log(Sig_j.determinant()) - std::log(Sig_i.determinant());
  return 0.5*result;
}

inline double nig_kl(const matrix& mu_i, const matrix& Sig_i, const real_type& nu_i, const real_type& iota_i,
                     const matrix& mu_j, const matrix& Sig_j, const real_type& nu_j, const real_type& iota_j)
{
  matrix I = matrix::Identity(Sig_i.rows(), Sig_i.cols());
  matrix Sig_j_inv = Sig_j.inverse();
  double result = (Sig_j_inv * Sig_i - I).trace();
  result = result + (nu_i/iota_i) * ((mu_j - mu_i).transpose() * Sig_j_inv * (mu_j - mu_i))(0,0);
  result = result + std::log(Sig_j.determinant()) - std::log(Sig_i.determinant());
  result = 0.5*result;
  result = result + nu_j * std::log(iota_i/iota_j);
  result = result - (std::lgamma(nu_i) - std::lgamma(nu_j));
  result = result + (nu_i - nu_j)*boost::math::digamma(nu_i);
  result = result - (iota_i - iota_j)*(nu_i/iota_i);
  return result;
}




template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<known_variance,KL_divergence_type::exact> >; }
{
  std::tuple<matrix,matrix> ta = transform(a);
  std::tuple<matrix,matrix> tb = transform(b);
  matrix mu_i  = std::get<0>(ta);
  matrix Sig_i = std::get<1>(ta);      
  matrix mu_j  = std::get<0>(tb);
  matrix Sig_j = std::get<1>(tb);
  real_type sigma = a.noise_structure.sigma;  

  real_type result;
  if(mu_i.rows() == 1)
    {
      real_type m_i = mu_i(0,0), m_j = mu_j(0,0);
      real_type V_i = sigma*sigma*Sig_i(0,0);   
      real_type V_j = sigma*sigma*Sig_j(0,0);
      if(std::abs(V_j - V_i) < 1e-12)
        {
          real_type s = std::sqrt(0.5*(V_i + V_j));
          result = 2.0*normal_cdf(std::abs(m_i - m_j)/(2.0*s)) - 1.0;
        }
      else
        {
          real_type aa = m_i*V_j - m_j*V_i;
          real_type bb = std::sqrt(V_i*V_j) *
                         std::sqrt((m_i - m_j)*(m_i - m_j) + (V_j - V_i)*std::log(V_j/V_i));
          real_type c1 = (aa - bb)/(V_j - V_i);
          real_type c2 = (aa + bb)/(V_j - V_i);
          if(c1 > c2){ real_type tmp = c1; c1 = c2; c2 = tmp; }
          result =  normal_cdf((c2 - m_i)/std::sqrt(V_i)) - normal_cdf((c1 - m_i)/std::sqrt(V_i))
                  + normal_cdf((c1 - m_j)/std::sqrt(V_j)) - normal_cdf((c2 - m_j)/std::sqrt(V_j));
        }
    }
  else
    {
      result = std::sqrt(0.5 * gaussian_kl(mu_i, Sig_i, mu_j, Sig_j, sigma));
    }
  return a.weight*result;

}

template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<known_variance,KL_divergence_type::approximate> >; }
{
  std::tuple<matrix,matrix> transformed = transform(a);
  matrix mu_i = std::get<0>(transformed);
  matrix sigma_i = std::get<1>(transformed);
  transformed = transform(b);
  matrix mu_i_plus_1 = std::get<0>(transformed);
  matrix sigma_i_plus_1 = std::get<1>(transformed);
  double sigma = a.noise_structure.sigma;
  return a.weight * std::sqrt(0.5 * gaussian_kl(mu_i, sigma_i, mu_i_plus_1, sigma_i_plus_1, sigma));
}

template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<unknown_variance,KL_divergence_type::exact> >; }
{
  std::tuple<matrix,matrix> transformed = transform(a);
  matrix mu_i = std::get<0>(transformed);
  matrix sigma_i = std::get<1>(transformed);
  transformed = transform(b);
  matrix mu_i_plus_1 = std::get<0>(transformed);
  matrix sigma_i_plus_1 = std::get<1>(transformed);
  return a.weight * std::sqrt(0.5 * nig_kl(mu_i, sigma_i, a.noise_structure.nu, a.noise_structure.iota,
                           mu_i_plus_1, sigma_i_plus_1, b.noise_structure.nu, b.noise_structure.iota));
}

template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<unknown_variance,KL_divergence_type::approximate> >; }
{
  std::tuple<matrix,matrix> transformed = transform(a);
  matrix mu_i = std::get<0>(transformed);
  matrix sigma_i = std::get<1>(transformed);
  transformed = transform(b);
  matrix mu_i_plus_1 = std::get<0>(transformed);
  matrix sigma_i_plus_1 = std::get<1>(transformed);
  return a.weight * std::sqrt(0.5 * nig_kl(mu_i, sigma_i, a.noise_structure.nu, a.noise_structure.iota,
                           mu_i_plus_1, sigma_i_plus_1, b.noise_structure.nu, b.noise_structure.iota));
}


template<typename noise, KL_divergence_type KL_divergence>
std::list<particle_type<noise,KL_divergence> >& prune(
  std::list<particle_type<noise,KL_divergence> >& particles,
  const int& max_num_particles,
  const pruning_type& pruning_strategy)
{
  if(particles.size() > max_num_particles && particles.size() > 3)
    {
      if(pruning_strategy == pruning_type::LP_P)
        { 
          auto it_first_candidate = particles.begin();
          it_first_candidate++;
          auto it_evicted = std::min_element(
            it_first_candidate,
            particles.end(),
            [](const auto& a, const auto& b)
            {
              return a.weight < b.weight;
            }
          );
          particles.erase(it_evicted);
          return particles;
        }
      auto it_1 = particles.begin();
      it_1++;
      auto it_2 = particles.begin();
      it_2++;it_2++;
      auto it_n_minus_1 = particles.end();
      it_n_minus_1--;
      std::list<double> total_variations;      
      std::transform(it_1,it_n_minus_1,it_2,
        std::back_inserter(total_variations),
        [&pruning_strategy](auto& a,auto& b){
          if (pruning_strategy == pruning_type::LP_M)
            {
              return static_cast<double>(a.weight);
            }
          return total_variation(a,b);});
      auto it_min_total_variation = std::min_element(total_variations.begin(),total_variations.end());
      auto it_evicted = particles.begin();
      std::advance(it_evicted,std::distance(total_variations.begin(),it_min_total_variation));
      auto it_relocation = it_evicted;
      it_evicted++; 
      it_relocation++;it_relocation++;

      auto combined_weight = it_relocation -> weight + it_evicted -> weight;
      if(combined_weight == 0.0)
      {
        std::fill(
          it_evicted->ratios.begin(),
          it_evicted->ratios.end(),
          0.0
        );

        std::fill(
          it_relocation->ratios.begin(),
          it_relocation->ratios.end(),
          0.0
        );
      }
      else
      {
        auto weight = it_evicted->weight;

        std::transform(it_evicted->ratios.begin(),
          it_evicted->ratios.end(),
          it_evicted->ratios.begin(),
          [&combined_weight,&weight](auto& ratio){return ratio*weight/combined_weight;});
        weight = it_relocation -> weight;
        std::transform(it_relocation->ratios.begin(),
          it_relocation->ratios.end(),
          it_relocation->ratios.begin(),
          [&combined_weight,&weight](auto& ratio){return ratio*weight/combined_weight;});}
      it_relocation->ratios.insert(it_relocation->ratios.begin(),it_evicted->ratios.begin(),it_evicted->ratios.end());
      it_relocation -> weight += it_evicted -> weight;
      particles.erase(it_evicted);
    }  
  return particles;
  
}



#endif
