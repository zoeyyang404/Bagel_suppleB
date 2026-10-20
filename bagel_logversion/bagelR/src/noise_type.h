
#ifndef ___NOISE_TYPE_H___
#define ___NOISE_TYPE_H___


struct known_variance 
{
  double sigma;
};


struct unknown_variance
{
  double nu;
  double iota;

  unknown_variance& operator=(const unknown_variance&);
};









#endif
