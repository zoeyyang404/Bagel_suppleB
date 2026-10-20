
setClass("bagel_type")


sequential_bagel_kv_exact <- function(feature_vector,prior,transform,p0,p,sigma,max_particles)
{
   bagel_object <- new("bagel_kv_exact_type",feature_vector=feature_vector,
				             prior=prior,
			                     transform=transform,
			                     p0=p0,
			                     p=p,
			                     sigma=sigma,
			                     max_particles=max_particles,
			                     bagel_object=new(bagelR_kv_exact,p0,p,sigma,max_particles))
  return(bagel_object)

}

sequential_bagel_uv_exact <- function(feature_vector,prior,transform,p0,p,nu,iota,max_particles)
{
   bagel_object <- new("bagel_uv_exact_type",feature_vector=feature_vector,
						   prior=prior,
			                  	   transform=transform,
			                  	   p0=p0,
			                  	   p=p,
			                  	   nu=nu,
			                  	   iota=iota,
			                  	   max_particles=max_particles,
			                  	   bagel_object=new(bagelR_uv_exact,p0,p,nu,iota,max_particles))
  return(bagel_object)
}


sequential_bagel_kv_approximate <- function(feature_vector,prior,transform,p0,p,sigma,max_particles)
{
   bagel_object <- new("bagel_kv_approximate_type",feature_vector=feature_vector,
				             prior=prior,
			                     transform=transform,
			                     p0=p0,
			                     p=p,
			                     sigma=sigma,
			                     max_particles=max_particles,
			                     bagel_object=new(bagelR_kv_approximate,p0,p,sigma,max_particles))
  return(bagel_object)

}

sequential_bagel_uv_approximate <- function(feature_vector,prior,transform,p0,p,nu,iota,max_particles)
{
   bagel_object <- new("bagel_uv_approximate_type",feature_vector=feature_vector,
						   prior=prior,
			                  	   transform=transform,
			                  	   p0=p0,
			                  	   p=p,
			                  	   nu=nu,
			                  	   iota=iota,
			                  	   max_particles=max_particles,
			                  	   bagel_object=new(bagelR_uv_approximate,p0,p,nu,iota,max_particles))
  return(bagel_object)				   
}




# set_weights
setGeneric("set_weights",function(object,weights) standardGeneric("set_weights"))
setMethod("set_weights","bagel_type",
          function(object,weights)
	  {
         object@bagel_object$set_weights(weights)     
	  })


# weights
setGeneric("weights",function(object) standardGeneric("weights"))
setMethod("weights","bagel_type",
          function(object)
	  {

	    return(Reduce(c,Map("*",object@bagel_object$get_weights(),object@bagel_object$get_ratios())))

	  })


# active_weights
setGeneric("active_weights",function(object) standardGeneric("active_weights"))
setMethod("active_weights","bagel_type",
          function(object)
	  {
        return(object@bagel_object$get_weights())
	  })

# active_posteriors
# The returned list has the same order and length as active_weights().
setGeneric("active_posteriors",function(object) standardGeneric("active_posteriors"))
setMethod("active_posteriors","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_posteriors())
	  })


# taus
setGeneric("taus",function(object) standardGeneric("taus"))
setMethod("taus","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_taus())
	  })



# ratios
setGeneric("ratios",function(object) standardGeneric("ratios"))
setMethod("ratios","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_ratios())	
	  })


# weights_tau_zero
setGeneric("weights_tau_zero",function(object) standardGeneric("weights_tau_zero"))
setMethod("weights_tau_zero","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_weights()[1])
	  })


# log_weights_tau_zero - no-change posterior weight on the log scale; stays
# finite and informative even when exp() of it would underflow to zero
setGeneric("log_weights_tau_zero",function(object) standardGeneric("log_weights_tau_zero"))
setMethod("log_weights_tau_zero","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_log_weights()[1])
	  })



# time
setGeneric("time",function(object) standardGeneric("time"))
setMethod("time","bagel_type",
          function(object)
	  {
	    # note - time is incremented preemptivley by bagel object
	    return(object@bagel_object$get_time() - 1)
	  })

# max_particles
setGeneric("max_particles",function(object) standardGeneric("max_particles"))
setMethod("max_particles","bagel_type",
          function(object)
	  {
	    return(object@bagel_object$get_max_num_particles())
	  })





setClass("bagel_result_type", slots=list(y = "numeric",
                                         threshold = "numeric",
					 max_particles = "numeric",
                                         w0ts = "numeric",
                                         log_w0ts = "numeric",
					 weights="numeric",
					 posteriors="list",
					 ratios="list",
					 taus = "numeric",
					 t="numeric"))

bagel_result <- function(y,
                         threshold,
			 max_particles,
                         w0ts,
			 log_w0ts,
			 weights,
			 posteriors,
			 ratios,
			 taus,
			 t)
{
	return(new("bagel_result_type",y=y,
	                               threshold=threshold,
				       max_particles=max_particles,
	                               w0ts=w0ts,
	                               log_w0ts=log_w0ts,
				       weights=weights,
				       posteriors=posteriors,
				       ratios=ratios,
				       taus=taus,
				       t=t))
}







setClass("bagel_uv_approximate_type", slots=list(feature_vector = "function",
      		       			         prior = "function",
				                 transform = "function",
				                 p0 = "numeric",
				                 p = "numeric",
				                 nu = "numeric",
	                                         iota = "numeric",	
				                 max_particles = "numeric",
				                 bagel_object = "ANY"),
						 contains="bagel_type")

bagel_uv_approximate <- function(Y,feature_vector,prior,transform,p0,p,nu,iota,max_particles,threshold,log_threshold=NA)
{

   bagel_object <- new("bagel_uv_approximate_type",feature_vector=feature_vector,
						   prior=prior,
			                  	   transform=transform,
			                  	   p0=p0,
			                  	   p=p,
			                  	   nu=nu,
			                  	   iota=iota,
			                  	   max_particles=max_particles,
			                  	   bagel_object=new(bagelR_uv_approximate,p0,p,nu,iota,max_particles))
  return(analyse(Y,threshold,bagel_object,log_threshold))

}

setClass("bagel_uv_exact_type", slots=list(feature_vector = "function",
      		       			   prior = "function",
				           transform = "function",
				           p0 = "numeric",
				           p = "numeric",
				           nu = "numeric",
	                                   iota = "numeric",	
				           max_particles = "numeric",
				           bagel_object = "ANY"),
					   contains="bagel_type")

bagel_uv_exact <- function(Y,feature_vector,prior,transform,p0,p,nu,iota,max_particles,threshold,log_threshold=NA)
{

   bagel_object <- new("bagel_uv_exact_type",feature_vector=feature_vector,
				             prior=prior,
			                     transform=transform,
			                     p0=p0,
			                     p=p,
			                     nu=nu,
			                     iota=iota,
			                     max_particles=max_particles,
			                     bagel_object=new(bagelR_uv_exact,p0,p,nu,iota,max_particles))
  return(analyse(Y,threshold,bagel_object,log_threshold))
}



setClass("bagel_kv_approximate_type", slots=list(feature_vector = "function",
      		       			         prior = "function",
				                 transform = "function",
				                 p0 = "numeric",
				                 p = "numeric",
				                 sigma = "numeric",
				                 max_particles = "numeric",
				                 bagel_object = "ANY"),
						 contains="bagel_type")

bagel_kv_approximate <- function(Y,feature_vector,prior,transform,p0,p,sigma,max_particles,threshold,log_threshold=NA)
{

   bagel_object <- new("bagel_kv_approximate_type",feature_vector=feature_vector,
						   prior=prior,
			                  	   transform=transform,
			                  	   p0=p0,
			                  	   p=p,
			                  	   sigma=sigma,
			                  	   max_particles=max_particles,
			                  	   bagel_object=new(bagelR_kv_approximate,p0,p,sigma,max_particles))
  return(analyse(Y,threshold,bagel_object,log_threshold))

}


setClass("bagel_kv_exact_type", slots=list(feature_vector = "function",
      		       			   prior = "function",
				           transform = "function",
				           p0 = "numeric",
				           p = "numeric",
				           sigma = "numeric",
				           max_particles = "numeric",
				           bagel_object = "ANY"),
					   contains="bagel_type")

bagel_kv_exact <- function(Y,feature_vector,prior,transform,p0,p,sigma,max_particles,threshold,log_threshold=NA)
{

   bagel_object <- new("bagel_kv_exact_type",feature_vector=feature_vector,
				             prior=prior,
			                     transform=transform,
			                     p0=p0,
			                     p=p,
			                     sigma=sigma,
			                     max_particles=max_particles,
			                     bagel_object=new(bagelR_kv_exact,p0,p,sigma,max_particles))
  return(analyse(Y,threshold,bagel_object,log_threshold))

}



# call

analyse <- function(Y,threshold,bagel_object,log_threshold=NA)
{
   # alarm on the log scale: 1 - w0 > threshold  <=>  -log(w0) > -log(1 - threshold).
   # The log form stays usable when w0 underflows to 0 (where 1 - w0 rounds to
   # exactly 1 and the probability-scale rule becomes uninformative).
   # threshold >= 1 (incl. Inf) is the established "never alarm" sentinel:
   # under the old rule 1 - w0 > threshold could not fire, so map it to Inf
   # instead of feeding log1p a value outside its domain.
   #
   # log_threshold: when supplied (not NA), it IS the alarm bound on the
   # -log(w0) scale and is used directly.  This bypasses the probability-scale
   # interface, whose 1 - exp(-x) rounds to 1 (=> never-alarm) for x > ~36.7,
   # so that thresholds beyond that ceiling (e.g. sigma2=1's ~50) can be applied.
   if (is.na(log_threshold))
     log_threshold <- if (threshold >= 1) Inf else -log1p(-threshold)
   log_w0ts <- list()
   for(y in Y)
   {
      update(bagel_object,y)
      log_w0t <- bagel_object@bagel_object$get_log_w0()
      log_w0ts <- append(log_w0ts,log_w0t)
      if(-log_w0t > log_threshold)
      {
	break
      }
   }
   log_w0ts <- as.numeric(log_w0ts)
   results <- bagel_result(as.numeric(Y),
                       threshold,
		       bagel_object@bagel_object$get_max_num_particles(),
                       exp(log_w0ts),
                       log_w0ts,
                       as.numeric(bagel_object@bagel_object$get_weights()),
		       bagel_object@bagel_object$get_posteriors(),
		       bagel_object@bagel_object$get_ratios(),
		       as.numeric(bagel_object@bagel_object$get_taus()),
                       bagel_object@bagel_object$get_time())
  return(results)

}












update <- function(object,y)
{
        bagel_object <- object@bagel_object 
	feature_vector <- object@feature_vector
	prior <- object@prior
	transform <- object@transform
	taus <- bagel_object$get_taus()
	t <- bagel_object$get_time()
	taus <- c(taus,t-1)
	bagel_object$set_feature_vectors(taus,Map(function(tau) return(feature_vector(t,tau)),taus))
	ts <- c(1,taus[-1],t)
	priors <- Map(prior,ts)
	prior_mus <- Map(function(x) return(x$mu), priors)
	prior_sigmas <- Map(function(x) x$sigma, priors)
	bagel_object$set_priors(ts,prior_mus,prior_sigmas)
	bagel_object$set_transformations(taus,Map(function(tau) return(transform(tau)),taus))
        result <- bagel_object$update(y)
	return(result)
}


# weights
#setGeneric("weights",function(object) standardGeneric("weights"))
setMethod("weights","bagel_result_type",
          function(object)
	  {
	    return(Reduce(c,Map("*",object@weights,object@ratios)))
	  })


# active_weights
#setGeneric("active_weights",function(object) standardGeneric("active_weights"))
setMethod("active_weights","bagel_result_type",
          function(object)
	  {
	    return(object@weights)
	  })

# active_posteriors
# The returned list has the same order and length as active_weights().
setMethod("active_posteriors","bagel_result_type",
          function(object)
	  {
	    return(object@posteriors)
	  })


# taus
#setGeneric("taus",function(object) standardGeneric("taus"))
setMethod("taus","bagel_result_type",
          function(object)
	  {
	    return(object@taus)
	  })



# ratios
#setGeneric("ratios",function(object) standardGeneric("ratios"))
setMethod("ratios","bagel_result_type",
          function(object)
	  {
	    return(object@ratios)	
	  })


# weights_tau_zero
#setGeneric("weights_tau_zero",function(object) standardGeneric("weights_tau_zero"))
setMethod("weights_tau_zero","bagel_result_type",
          function(object)
	  {
	    return(object@w0ts)
	  })


# log_weights_tau_zero
#setGeneric("log_weights_tau_zero",function(object) standardGeneric("log_weights_tau_zero"))
setMethod("log_weights_tau_zero","bagel_result_type",
          function(object)
	  {
	    return(object@log_w0ts)
	  })



# time
#setGeneric("time",function(object) standardGeneric("time"))
setMethod("time","bagel_result_type",
          function(object)
	  {
	    # note - time is incremented preemptivley by bagel object
	    return(object@t - 1)
	  })

# max_particles
#setGeneric("max_particles",function(object) standardGeneric("max_particles"))
setMethod("max_particles","bagel_result_type",
          function(object)
	  {
	    return(object@max_particles)
	  })


# changepoint
setGeneric("changepoint",function(object) standardGeneric("changepoint"))
setMethod("changepoint","bagel_result_type",
          function(object)
	  {
	    if(object@t - 1 < length(object@y))
	    {
	      return(list("location"=which.max(weights(object)[-1]),"detected"=time(object)))
	    }
	    return(NA)
	  })



