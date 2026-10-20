
#include "prior_type.h"


prior_type& prior_type::operator=(const prior_type& other)
{
  this -> mu = other.mu;
  this -> sigma = other.sigma;
  return *this;
}

