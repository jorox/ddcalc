#ifndef CAI_H
#define CAI_H

#include "segment.h"

namespace cai{
  double self_energy (double mu, double nu, double cw,
                      const Segment& neo);
  double interaction_energy ( double mu, double nu, double cw,
                              const Segment& neo, const Segment& other);
};

#endif
