#ifndef NODE_H
#define NODE_H

#include "Eigen/Dense"
#include <vector>

//class Segment2;

class Node
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Node( Eigen::Vector3d& coords );
  ~Node();
  //Node * connect (Segment2* seg);
  //Node * disconnect (Segment2* seg);
  Eigen::Vector3d get_coords() const;
  void change_coords(Eigen::Vector3d& newCoords);

 private:
  Eigen::Vector3d _x;
  //std::vector<Segment2*> _connectedSegments;
};
#endif
