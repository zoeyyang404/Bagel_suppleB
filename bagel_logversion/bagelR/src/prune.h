
#ifndef ___PRUNE_H___
#define ___PRUNE_H___


#include "particle_type.h"
#include "KL_divergence_type.h"
#include "normal_cdf.h"
#include <list>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <limits>
#include <algorithm>

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

// log(exp(a) + exp(b)) computed without leaving the log scale
inline real_type log_add_exp(const real_type& a,const real_type& b)
{
  if(a == -std::numeric_limits<real_type>::infinity()) return b;
  if(b == -std::numeric_limits<real_type>::infinity()) return a;
  real_type m = std::max(a,b);
  return m + std::log1p(std::exp(std::min(a,b) - m));
}

// KL divergence between two Gaussians N(mu_i, sigma^2 Sig_i) and N(mu_j, sigma^2 Sig_j)
// (Theorem 5). Shared by the known-variance 'approximate' branch and the
// multivariate fall-back of the 'exact' branch so they are bit-for-bit identical.
inline double gaussian_kl(const matrix& mu_i, const matrix& Sig_i,
                          const matrix& mu_j, const matrix& Sig_j,
                          const real_type& sigma)
{
  matrix I = matrix::Identity(Sig_i.rows(), Sig_i.cols());
  matrix Sig_j_inv = Sig_j.inverse();          // once; was recomputed on the next line
  double result = (Sig_j_inv * Sig_i - I).trace();
  result = result + (1.0/(sigma*sigma)) * ((mu_j - mu_i).transpose() * Sig_j_inv * (mu_j - mu_i))(0,0);
  // log det(Sig_i^-1 Sig_j) == log det(Sig_j) - log det(Sig_i), since
  // det(AB) = det(A) det(B) and det(A^-1) = 1/det(A). Computing it this way
  // avoids a matrix inversion and a matrix-matrix product.
  result = result + std::log(Sig_j.determinant()) - std::log(Sig_i.determinant());
  return 0.5*result;
}

// KL divergence between two Normal-Inverse-Gamma distributions (Theorem 6).
// Used for the unknown-variance case for BOTH the 'exact' and 'approximate'
// selectors -- there is no closed-form total variation for NIG, so unknown
// variance always uses this KL bound.
inline double nig_kl(const matrix& mu_i, const matrix& Sig_i, const real_type& nu_i, const real_type& iota_i,
                     const matrix& mu_j, const matrix& Sig_j, const real_type& nu_j, const real_type& iota_j)
{
  matrix I = matrix::Identity(Sig_i.rows(), Sig_i.cols());
  matrix Sig_j_inv = Sig_j.inverse();          // once; was recomputed on the next line
  double result = (Sig_j_inv * Sig_i - I).trace();
  result = result + (nu_i/iota_i) * ((mu_j - mu_i).transpose() * Sig_j_inv * (mu_j - mu_i))(0,0);
  // -log(det(Sig_i)/det(Sig_j)) == log det(Sig_j) - log det(Sig_i); same value,
  // written as a difference of logs (matches gaussian_kl and is more robust).
  result = result + std::log(Sig_j.determinant()) - std::log(Sig_i.determinant());
  result = 0.5*result;
  result = result + nu_j * std::log(iota_i/iota_j);
  result = result - (std::lgamma(nu_i) - std::lgamma(nu_j));
  result = result + (nu_i - nu_j)*boost::math::digamma(nu_i);
  result = result - (iota_i - iota_j)*(nu_i/iota_i);
  return result;
}




// ******************************************************************************
template <typename T>
double total_variation(const T& a,const T& b)
requires requires { requires std::same_as<T,particle_type<known_variance,KL_divergence_type::exact> >; }
{
  // Merge criterion for the post-change parameter theta (known noise variance):
  //   * theta is 1-D  -> EXACT total variation (Proposition 1 / Corollary 1)
  //   * theta is >1-D -> no closed-form TV exists, so fall back to the KL
  //                      (Pinsker) bound (Theorem 5), as in the 'approximate' branch.
  // The dimension of theta is the number of rows of the transformed mean.
  std::tuple<matrix,matrix> ta = transform(a);
  std::tuple<matrix,matrix> tb = transform(b);
  matrix mu_i  = std::get<0>(ta);
  matrix Sig_i = std::get<1>(ta);          // scaled covariance of theta
  matrix mu_j  = std::get<0>(tb);
  matrix Sig_j = std::get<1>(tb);
  real_type sigma = a.noise_structure.sigma;   // known noise sd

  real_type result;
  if(mu_i.rows() == 1)
    {
      // ----- exact 1-D total variation -----
      real_type m_i = mu_i(0,0), m_j = mu_j(0,0);
      real_type V_i = sigma*sigma*Sig_i(0,0);   // actual variance = noise var * scaled cov
      real_type V_j = sigma*sigma*Sig_j(0,0);
      if(std::abs(V_j - V_i) < 1e-12)
        {
          // equal-variance closed form: TV = 2 Phi(|dmu|/(2 sqrt(V))) - 1
          real_type s = std::sqrt(0.5*(V_i + V_j));
          result = 2.0*normal_cdf(std::abs(m_i - m_j)/(2.0*s)) - 1.0;
        }
      else
        {
          // the two densities intersect at c1 < c2
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
      // ----- multivariate fall-back: KL divergence (Theorem 5), shared helper -----
      result = std::sqrt(0.5 * gaussian_kl(mu_i, Sig_i, mu_j, Sig_j, sigma));
    }
  return std::exp(a.log_weight)*result;

}
// ******************************************************************************

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
  return std::exp(a.log_weight) * std::sqrt(0.5 * gaussian_kl(mu_i, sigma_i, mu_i_plus_1, sigma_i_plus_1, sigma));
}

// ******************************************************************************
// Unknown variance: NIG KL bound (Theorem 6) for BOTH selectors. There is no
// closed-form total variation for the Normal-Inverse-Gamma posterior, so the
// 'exact' and 'approximate' branches are identical (both use the KL bound).
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
  return std::exp(a.log_weight) * std::sqrt(0.5 * nig_kl(mu_i, sigma_i, a.noise_structure.nu, a.noise_structure.iota,
                           mu_i_plus_1, sigma_i_plus_1, b.noise_structure.nu, b.noise_structure.iota));
}
// ******************************************************************************


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
  return std::exp(a.log_weight) * std::sqrt(0.5 * nig_kl(mu_i, sigma_i, a.noise_structure.nu, a.noise_structure.iota,
                           mu_i_plus_1, sigma_i_plus_1, b.noise_structure.nu, b.noise_structure.iota));
}


template<typename noise, KL_divergence_type KL_divergence>
std::list<particle_type<noise,KL_divergence> >& prune(std::list<particle_type<noise,KL_divergence> >& particles,const int& max_num_particles)
{
  // return particles;
  if(particles.size() > max_num_particles && particles.size() > 3)
    {
      // auto it_1 = particles.begin();
      // it_1++;
      // auto it_2 = particles.begin();
      // it_2++;it_2++;
      // auto it_n_minus_1 = particles.end();
      // it_n_minus_1--;
      // // generate total variations between adjacent particles
      // std::list<double> total_variations;      
      // std::transform(it_1,it_n_minus_1,it_2,std::back_inserter(total_variations),[](auto& a,auto& b){return total_variation(a,b);});
      // // locate the minimum total variation
      // auto it_min_total_variation = std::min_element(total_variations.begin(),total_variations.end());
      // // locate the particle to be evicted from the particle population
      // auto it_evicted = particles.begin();
      // std::advance(it_evicted,std::distance(total_variations.begin(),it_min_total_variation));
      // auto it_relocation = it_evicted;
      // it_evicted++; 
      // it_relocation++;it_relocation++;
      auto it_first_cp = particles.begin();
      it_first_cp++;  // skip tau = 0

      auto it_second_cp = it_first_cp;
      it_second_cp++;

      auto it_newest = particles.end();
      it_newest--;    // the candidate introduced at this iteration (tau = t-1)

      // Compare EVERY adjacent pair of changepoint particles, including the
      // pair whose right-hand element is the newly introduced candidate:
      //   (first_cp, second_cp), ..., (second_newest, newest)
      // std::transform pairs [it_first_cp, it_newest) with [it_second_cp, ...),
      // so the LEFT element ranges over every changepoint particle except the
      // newest.  The newest can therefore ABSORB its predecessor but is never
      // itself evicted -- prune always removes the left element of the chosen
      // pair and relocates its weight to the right neighbour.  The tau = 0
      // particle is never a merge candidate.
      std::list<double> total_variations;

      std::transform(
        it_first_cp,
        it_newest,
        it_second_cp,
        std::back_inserter(total_variations),
        [](auto& a, auto& b) {
          return total_variation(a, b);
        }
      );

      if (total_variations.empty()) {
        return particles;
      }

      // locate the minimum total variation
      auto it_min_total_variation =
        std::min_element(total_variations.begin(), total_variations.end());

      // locate the particle to be evicted
      auto it_evicted = it_first_cp;
      std::advance(
        it_evicted,
        std::distance(total_variations.begin(), it_min_total_variation)
      );

      // relocate weight to the right neighbour
      auto it_relocation = it_evicted;
      it_relocation++;
      // combine the two weights on the log scale; the ratio rescaling factors
      // weight/combined_weight are formed as exp of log differences, which are
      // well conditioned even when both weights have underflowed linear scale
      auto log_combined_weight = log_add_exp(it_relocation -> log_weight,it_evicted -> log_weight);
      auto weight_fraction = std::exp(it_evicted -> log_weight - log_combined_weight);
      std::transform(it_evicted->ratios.begin(),
		     it_evicted->ratios.end(),
		     it_evicted->ratios.begin(),
		     [&weight_fraction](auto& ratio){return ratio*weight_fraction;});
      weight_fraction = std::exp(it_relocation -> log_weight - log_combined_weight);
      std::transform(it_relocation->ratios.begin(),
		     it_relocation->ratios.end(),
		     it_relocation->ratios.begin(),
		     [&weight_fraction](auto& ratio){return ratio*weight_fraction;});
      it_relocation->ratios.insert(it_relocation->ratios.begin(),it_evicted->ratios.begin(),it_evicted->ratios.end());
      // update weight
      it_relocation -> log_weight = log_combined_weight;
      // evict the pruned particle
      particles.erase(it_evicted);
    }  
  return particles;
  
}



#endif
