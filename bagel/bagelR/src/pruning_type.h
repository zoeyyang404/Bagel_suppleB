#ifndef BAGEL_PRUNING_TYPE_H
#define BAGEL_PRUNING_TYPE_H

#include <stdexcept>
#include <string>

enum class pruning_type
{
  bagel,
  LP_M,
  LP_P
};

inline pruning_type parse_pruning_type(const std::string& value)
{
  if(value == "bagel")
    {
      return pruning_type::bagel;
    }

  if(value == "LP-M")
    {
      return pruning_type::LP_M;
    }

  if(value == "LP-P")
    {
      return pruning_type::LP_P;
    }

  throw std::invalid_argument(
    "pruning must be one of: 'bagel', 'LP-M', or 'LP-P'"
  );
}

#endif