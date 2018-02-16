#ifndef SEGMENT_H
#define SEGMENT_H

#include "Eigen/Dense"

class Segment{

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Segment (const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& brg);
  Segment (const Segment& );
  ~Segment();

  double length() const;
  void get_unit_vector(Eigen::Vector3d& ) const;
  void get_burgers_vector(Eigen::Vector3d& ) const;
  void get_head(Eigen::Vector3d& ) const;
  void get_tail(Eigen::Vector3d& ) const;
  //Segment& operator= (Segment const &rhs);

 private:
  Eigen::Vector3d _x1;
  Eigen::Vector3d _x2;
  Eigen::Vector3d _burgers;
  Eigen::Vector3d _unitVector;
  Eigen::Vector3d _lineVector;
};

#endif
