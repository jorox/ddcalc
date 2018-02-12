#ifndef CAI_H
#define CAI_H

#include "segment.h"

class Cai{
 public:
  Cai (double mu, double nu, double cw);
  double self_energy (const Segment& neo);
  double interaction_energy (const Segment& neo, const Segment& other);
  void calculate_total_energy(const std::vector<Segment>&, double*);

 private:
  double _mu;
  double _nu;
  double _cw;

 protected:
  double special_func ( const Eigen::Vector3d& R,
                        const Eigen::Vector3d& t1,
                        const Eigen::Vector3d& t2,
                        const Eigen::Vector3d& b1,
                        const Eigen::Vector3d& b2 ) const;

  double special_func_parallel( const Eigen::Vector3d& R,
                                const Eigen::Vector3d& t,
                                const Eigen::Vector3d& b1,
                                const Eigen::Vector3d& b2) const;
};

#endif
