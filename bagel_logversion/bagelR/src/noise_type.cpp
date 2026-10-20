#include "noise_type.h"


unknown_variance& unknown_variance::operator=(const unknown_variance& other)
  {
    nu = other.nu;
    iota = other.iota;
    return *this;
  }
