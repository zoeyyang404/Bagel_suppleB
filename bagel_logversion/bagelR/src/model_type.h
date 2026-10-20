#ifndef ___MODEL_TYPE_H___
#define ___MODEL_TYPE_H___

#include <functional>
#include "matrix_type.h"

#include "prior_type.h"
#include "time_type.h"


typedef matrix feature_vector_type;
typedef std::function<feature_vector_type (const int&,const int&)> feature_vector_function_type ;
typedef std::function<matrix (const time_type&)> transformation_function_type;
typedef std::function<prior_type (const int&)> prior_function_type;

struct model_type
{
  prior_function_type prior_function;
  feature_vector_function_type feature_vector_function;
  transformation_function_type transformer_function;
};


#endif
