#ifndef HELICALTURN_H
#define HELICALTURN_H

#include "segment.h"
#include "Eigen/Dense"
#include <stdio.h>
#include <vector>

class HelicalTurn
{
 public:
  HelicalTurn (double lmin, double totalLength, double height, double width, double burgers);
  //~HelicalTurn();
  void regenerate (double prismaticLength, double inclination);
  double return_inclination_angle() const;
  void get_segments( const std::vector<Segment>& ) const;
  void write_to_file(FILE * pf) const;

 private:
  double _angle;
  double _pLength;
  double _totalLength;
  double _height;
  double _width;
  Eigen::Vector3d _burgers;
  double _dx;
  std::vector<Segment> _segVector;

  void subdivide(const Segment& s1);
};
#endif
