
#ifndef ___PRIOR_TYPE_H___
#define ___PRIOR_TYPE_H___

#include <functional>
#include "matrix_type.h"
 

struct prior_type
{
  matrix mu;
  matrix sigma;
  prior_type& operator=(const prior_type&);
};


#endif
