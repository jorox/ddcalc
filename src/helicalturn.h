#ifndef HELICALTURN_H
#define HELICALTURN_H

#include "segment.h"
#include "Eigen/Dense"
#include <stdio.h>
#include <vector>

class HelicalTurn
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  HelicalTurn (double lmin, double totalLength, double height, double width, double burgers);
  //~HelicalTurn();
  void generate (double prismaticLength, double inclination, int image);
  //void get_segments( std::vector<Segment>& ) const;
  double calculate_minimum_angle(double prismaticLength) const;
  double get_total_length() const;
  void write_to_file(FILE * pf) const;
  std::vector<Segment> _segVector;

 private:
  double _angle;
  double _totalLength;
  double _height;
  double _width;
  Eigen::Vector3d _burgers;
  double _dx;

  void subdivide(const Segment& s1);
};
#endif
