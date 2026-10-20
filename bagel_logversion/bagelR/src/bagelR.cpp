
#include <Rcpp.h>
using namespace Rcpp;


#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

#include <algorithm>
#include <cmath>
#include <list>
#include <vector>
#include <map>
#include <memory>

#include "bagel.h"


template<typename noise, KL_divergence_type KL_divergence>
struct bagelR
{   
  std::shared_ptr<bagel_type<noise,KL_divergence> > sp_bagel; 
  std::map<int,matrix> M;
  std::map<int,prior_type> P;
  std::map<int,matrix> A;
  
  noise noise_structure;

  bagelR(const probability_type& p0,
	 const probability_type& p,
	 const real_type& nu,
	 const real_type& iota,
	 const int& n)
  {
    noise_structure.nu = nu;
    noise_structure.iota = iota;
    
    model_type model;
  
    model.prior_function =  std::bind(&bagelR<noise,KL_divergence>::prior_from_R, this, std::placeholders::_1);
    model.feature_vector_function =  std::bind(&bagelR<noise,KL_divergence>::feature_vector_from_R, this, std::placeholders::_1,std::placeholders::_2);
    // the transformer needs changing after methods have been added
    model.transformer_function =  std::bind(&bagelR<noise,KL_divergence>::transformation_from_R, this, std::placeholders::_1);
  
    sp_bagel = std::make_shared<bagel_type<noise,KL_divergence> >(bagel_type<noise,KL_divergence>(model,noise_structure,p0,p,n));  
  }
  
  bagelR(const probability_type& p0,
	 const probability_type& p,
	 const real_type& sigma,
	 const int& n)
  {
    noise_structure.sigma = sigma;
    
    model_type model;
  
    model.prior_function =  std::bind(&bagelR<noise,KL_divergence>::prior_from_R, this, std::placeholders::_1);
    model.feature_vector_function =  std::bind(&bagelR<noise,KL_divergence>::feature_vector_from_R, this, std::placeholders::_1,std::placeholders::_2);
    // the transformer needs changing after methods have been added
    model.transformer_function =  std::bind(&bagelR<noise,KL_divergence>::transformation_from_R, this, std::placeholders::_1);
    // model.transformer_function =  std::bind(&bagelR<noise,KL_divergence>::feature_vector_from_R, this, std::placeholders::_1,std::placeholders::_2);
  
    sp_bagel = std::make_shared<bagel_type<noise,KL_divergence> >(bagel_type<noise,KL_divergence>(model,noise_structure,p0,p,n));  
  }
  

  ~bagelR()
  {
  }
  
  real_type update(const real_type& x)
  {
    sp_bagel -> update(x);
    return sp_bagel -> weight_0_t();
  }

  // no-change posterior weight on the log scale, taken directly from the
  // internally retained log weight (no exponentiation, no precision loss)
  real_type get_log_w0()
  {
    return sp_bagel -> log_weight_0_t();
  }



void set_weights(const std::vector<real_type>& weights)
  {
    auto it_particle = sp_bagel->particles.begin();
    for(auto& weight : weights)
      {
	it_particle -> log_weight = std::log(weight);
	it_particle++;
      }
  }


  void set_log_weights(const std::vector<real_type>& log_weights)
  {
    auto it_particle = sp_bagel->particles.begin();
    for(auto& log_weight : log_weights)
      {
	it_particle -> log_weight = log_weight;
	it_particle++;
      }
  }


  std::list<double> get_weights()
  {

    std::list<double> lweights;
    std::transform(sp_bagel->particles.begin(),
		   sp_bagel->particles.end(),
		   std::inserter(lweights,lweights.end()),
		   [](auto& particle){return std::exp(particle.log_weight);});
    return lweights;
  }

  Rcpp::List get_posteriors()
  {
    Rcpp::List posteriors;
    for(const auto& particle : sp_bagel->particles)
      {
        Rcpp::List posterior = Rcpp::List::create(
          Rcpp::Named("tau") = particle.tau,
          Rcpp::Named("weight") = std::exp(particle.log_weight),
          Rcpp::Named("mu") = particle.post.mu,
          Rcpp::Named("sigma") = particle.post.sigma
        );
        if constexpr(std::same_as<noise,unknown_variance>)
          {
            posterior["nu"] = particle.noise_structure.nu;
            posterior["iota"] = particle.noise_structure.iota;
          }
        posteriors.push_back(posterior);
      }
    return posteriors;
  }


  std::list<double> get_log_weights()
  {

    std::list<double> llog_weights;
    std::transform(sp_bagel->particles.begin(),
		   sp_bagel->particles.end(),
		   std::inserter(llog_weights,llog_weights.end()),
		   [](auto& particle){return particle.log_weight;});
    return llog_weights;
  }
  
  double get_time()
  {
    return sp_bagel -> t;
  }
  
  
  std::list<int> get_taus()
  {
    std::list<int> ltaus;
    std::transform(sp_bagel->particles.begin(),
		   sp_bagel->particles.end(),
		   std::inserter(ltaus,ltaus.end()),
		   [](auto& particle){return particle.tau;});  
    return ltaus;  
  }
  
  const matrix& feature_vector_from_R(const int& t,const int& tau)
  { 
    return M[tau];
  }
  
  const prior_type prior_from_R(const int& t)
  {
    return P[t];
  }

  const matrix& transformation_from_R(const int& tau)
  {    
    return A[tau];
  }


  void set_transformations(const std::vector<int>& taus_from_R, const std::list<matrix>& transformations_from_R)
  {
    A.clear();
    std::transform(taus_from_R.begin(),
		   taus_from_R.end(),
		   transformations_from_R.begin(),
		   std::inserter(A,A.end()),
		   [](auto& tau,auto& a){return std::make_pair(tau,a);});
  }
  
  void set_feature_vectors(const std::vector<int>& taus_from_R, const std::list<matrix>& feature_vectors_from_R)
  { 
    M.clear();
    std::transform(taus_from_R.begin(),
		   taus_from_R.end(),
		   feature_vectors_from_R.begin(),
		   std::inserter(M,M.end()),
		   [](auto& tau,auto& m){return std::make_pair(tau,m);});
  }
  
  void set_priors(const std::vector<int>& ts_from_R,
		  const std::list<matrix>& prior_mus_from_R,
		  const std::list<matrix>& prior_sigmas_from_R)
  {
    P.clear();
    
    // zip into prior_types
    std::list<prior_type> priors;
    std::transform(prior_mus_from_R.begin(),
		   prior_mus_from_R.end(),
		   prior_sigmas_from_R.begin(),
		   std::inserter(priors,priors.end()),
		   [](auto& mu,auto& sigma){prior_type prior;
		     prior.mu = mu;
		     prior.sigma = sigma;
		     return prior;});
    std::transform(ts_from_R.begin(),
		   ts_from_R.end(),
		   priors.begin(),
		   std::inserter(P,P.end()),
		   [](auto& t,auto& prior){return std::make_pair(t,prior);});
  }
  
  
  
  std::list<std::list<double> > get_ratios()
  {
    return sp_bagel -> ratios();
  }

int get_max_num_particles()
  {
    return sp_bagel -> max_num_particles;
  }
  
  
};

typedef bagelR<unknown_variance,KL_divergence_type::exact> bagelR_uv_exact;
typedef bagelR<known_variance,KL_divergence_type::exact> bagelR_kv_exact;
typedef bagelR<unknown_variance,KL_divergence_type::approximate> bagelR_uv_approximate;
typedef bagelR<known_variance,KL_divergence_type::approximate> bagelR_kv_approximate;

RCPP_MODULE(bagelR) 
{
  class_<bagelR_uv_exact>("bagelR_uv_exact")
  .constructor<probability_type,probability_type,real_type,real_type,int>()
  .method("get_time", &bagelR_uv_exact::get_time)
  .method("get_taus", &bagelR_uv_exact::get_taus)
  .method("get_weights", &bagelR_uv_exact::get_weights)
  .method("get_posteriors", &bagelR_uv_exact::get_posteriors)
  .method("set_weights", &bagelR_uv_exact::set_weights)
  .method("get_log_weights", &bagelR_uv_exact::get_log_weights)
  .method("set_log_weights", &bagelR_uv_exact::set_log_weights)
  .method("get_log_w0", &bagelR_uv_exact::get_log_w0)
  .method("set_feature_vectors", &bagelR_uv_exact::set_feature_vectors)
  .method("set_priors", &bagelR_uv_exact::set_priors)
  .method("set_transformations", &bagelR_uv_exact::set_transformations)
  .method("update", &bagelR_uv_exact::update)
  .method("get_ratios", &bagelR_uv_exact::get_ratios)
  .method("get_max_num_particles", &bagelR_uv_exact::get_max_num_particles)
;
  class_<bagelR_uv_approximate>("bagelR_uv_approximate")
  .constructor<probability_type,probability_type,real_type,real_type,int>()
  .method("get_time", &bagelR_uv_approximate::get_time)
  .method("get_taus", &bagelR_uv_approximate::get_taus)
  .method("get_weights", &bagelR_uv_approximate::get_weights)
  .method("get_posteriors", &bagelR_uv_approximate::get_posteriors)
  .method("set_weights", &bagelR_uv_approximate::set_weights)
  .method("get_log_weights", &bagelR_uv_approximate::get_log_weights)
  .method("set_log_weights", &bagelR_uv_approximate::set_log_weights)
  .method("get_log_w0", &bagelR_uv_approximate::get_log_w0)
  .method("set_feature_vectors", &bagelR_uv_approximate::set_feature_vectors)
  .method("set_priors", &bagelR_uv_approximate::set_priors)
  .method("set_transformations", &bagelR_uv_approximate::set_transformations)
  .method("update", &bagelR_uv_approximate::update)
  .method("get_ratios", &bagelR_uv_approximate::get_ratios)
  .method("get_max_num_particles", &bagelR_uv_approximate::get_max_num_particles)
;



  class_<bagelR_kv_exact>("bagelR_kv_exact")
  .constructor<probability_type,probability_type,real_type,int>()
  .method("get_time", &bagelR_kv_exact::get_time)
  .method("get_taus", &bagelR_kv_exact::get_taus)
  .method("get_weights", &bagelR_kv_exact::get_weights)
  .method("get_posteriors", &bagelR_kv_exact::get_posteriors)
  .method("set_weights", &bagelR_kv_exact::set_weights)
  .method("get_log_weights", &bagelR_kv_exact::get_log_weights)
  .method("set_log_weights", &bagelR_kv_exact::set_log_weights)
  .method("get_log_w0", &bagelR_kv_exact::get_log_w0)
  .method("set_feature_vectors", &bagelR_kv_exact::set_feature_vectors)
  .method("set_priors", &bagelR_kv_exact::set_priors)
  .method("set_transformations", &bagelR_kv_exact::set_transformations)
  .method("update", &bagelR_kv_exact::update)
  .method("get_ratios", &bagelR_kv_exact::get_ratios)
  .method("get_max_num_particles", &bagelR_kv_exact::get_max_num_particles)
;
  class_<bagelR_kv_approximate>("bagelR_kv_approximate")
  .constructor<probability_type,probability_type,real_type,int>()
  .method("get_time", &bagelR_kv_approximate::get_time)
  .method("get_taus", &bagelR_kv_approximate::get_taus)
  .method("get_weights", &bagelR_kv_approximate::get_weights)
  .method("get_posteriors", &bagelR_kv_approximate::get_posteriors)
  .method("set_weights", &bagelR_kv_approximate::set_weights)
  .method("get_log_weights", &bagelR_kv_approximate::get_log_weights)
  .method("set_log_weights", &bagelR_kv_approximate::set_log_weights)
  .method("get_log_w0", &bagelR_kv_approximate::get_log_w0)
  .method("set_feature_vectors", &bagelR_kv_approximate::set_feature_vectors)
  .method("set_priors", &bagelR_kv_approximate::set_priors)
  .method("set_transformations", &bagelR_kv_approximate::set_transformations)
  .method("update", &bagelR_kv_approximate::update)
  .method("get_ratios", &bagelR_kv_approximate::get_ratios)
  .method("get_max_num_particles", &bagelR_kv_approximate::get_max_num_particles)
;

}



