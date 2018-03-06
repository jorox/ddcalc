#ifndef SEGMENT2_H
#define SEGMENT2_H

#include "node.h"
#include "Eigen/Dense"

class Segment2
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Segment2 (Node* p1, Node* p2, const Eigen::Vector3d& brg);
  //Segment2 (const Segment& );
  ~Segment2();

  double length();
  void get_unit_vector(Eigen::Vector3d& );
  void get_burgers_vector(Eigen::Vector3d& );
  Node* head_node() const;
  Node* tail_node() const;

 private:
  Node* _n1;
  Node* _n2;
  Eigen::Vector3d _x1;
  Eigen::Vector3d _x2;
  Eigen::Vector3d _burgers;
  Eigen::Vector3d _unitVector;
  Eigen::Vector3d _lineVector;
  void refresh_data();
};

#endif
