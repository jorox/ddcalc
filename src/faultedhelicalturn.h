#ifndef FAULTEDHELICALTURN_H
#define FAULTEDHELICALTURN_H

#include "node.h"
#include "segment2.h"
#include "Eigen/Dense"
#include <stdio.h>
#include <vector>

class FaultedHelicalTurn
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  FaultedHelicalTurn (double lmin, double totalLength, double height, double width, double burgers);
  //~HelicalTurn();
  void generate (double prismaticLength, double inclination, int image);
  double calculate_minimum_angle(double prismaticLength) const;
  double get_total_length();
  //void write_to_file(FILE * pf) const;
  void write_to_xml(FILE * pf);
  int find_node_tag(Node* targetNode);
  std::vector<Segment2> _segVector;
  std::vector<Node> _nodeVector1;
  std::vector<Node> _nodeVector2;

 private:
  double _angle;
  double _totalLength;
  double _height;
  double _width;
  Eigen::Vector4i _burgersp;
  Eigen::Vector4i _burgersb1;
  Eigen::Vector4i _burgersb2;

  double _dx;

};
#endif
