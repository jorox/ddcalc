#include "helicalturn.h"
#include "segment.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#define _USE_MATH_DEFINES

HelicalTurn::HelicalTurn( double lmin,
                          double totalLength,
                          double height,
                          double width, double b)
{
  //
  this->_dx = lmin;
  this->_height = height;
  this->_width = width;
  this->_totalLength = totalLength;
  this->_burgers = b * Eigen::Vector3d(1.0,0.0,0.0);

  this->regenerate(this->_totalLength / 2.0, M_PI / 2.0, 0);
}

void HelicalTurn::subdivide(const Segment &s1)
{
  double dx = this->_dx;
  int n = (int) (s1.length() / dx);
  double rem = (s1.length() - (double) (n * dx));
  Eigen::Vector3d start;
  Eigen::Vector3d end;
  Eigen::Vector3d t;
  s1.get_head(start);
  s1.get_tail(end);
  s1.get_unit_vector(t);


  // psuh first segment
  this->_segVector.push_back(Segment(start, start + t * dx, this->_burgers));
  const Segment * p0 = &(this->_segVector.back());
  //printf( "DEBUG: start = %f %f %f\n", start(0), start (1), start(2));

  // push middle segments
  for (int i = 1; i < n-1; ++i){
    this->_segVector.back().get_tail(start);
    this->_segVector.push_back(Segment(start, start + t * dx, this->_burgers));
    //printf( "DEBUG: start = %f %f %f\n", start(0), start (1), start(2));
  }

  // push last segment which will probably have length < dx
  this->_segVector.back().get_tail(start);
  //printf( "DEBUG: start = %f %f %f\n", start(0), start (1), start(2));
  this->_segVector.push_back(Segment(start, start + rem * t, this->_burgers));

  p0->get_head(start);
  this->_segVector.back().get_tail(end);

  printf ("   generated %i segments, (%f, %f, %f) --> (%f, %f, %f)\n",
          n, start(0), start(1), start(2), end(0), end(1), end(2));
}

void HelicalTurn::regenerate(double pLength,
                             double inclination,
                             int image = 0)
{
  double minInc = this->calculate_minimum_angle(pLength);
  printf( "prismatic segment = %f, min angle = %f\n", pLength, minInc * 180. / M_PI);
  if (inclination < minInc) {
    printf ( "ERROR: HelicalTurn::regenerate, minimum possible angle = %f degrees\n", minInc * 180. / M_PI);
    return;
  }

  _segVector.clear();

  double tn = tan( inclination );
  double ho2 = this->_height / 2.0;
  double ho2otn = this->_height / 2.0 / tn;
  double lmini = (this->_totalLength - pLength - 4.0 * ho2otn) / 2.0;
  printf ("DEBUG: lmin = %f\n", lmini);


  // generate first prismatic segment
  int nMainNodes = 8 + 7 * 2 * image;
  Eigen::Matrix<double, 3, 8> x;
  Eigen::Matrix<double, 3, 1> pbcShift;
  int ix = 0;

  for (int ipbc = -image; ipbc < image + 1; ++ipbc){
    ix = ipbc + image;
    pbcShift << double(ipbc) * this->_totalLength, 0.0, 0.0;

    x.block<3,1>(0,ix) << -this->_totalLength / 2.0 , -this->_width, 0.0;
    x.block<3,1>(0,ix+1) << x.block<3,1>(0,0) + Eigen::Vector3d(pLength/2.0, 0.0, 0.0);
    x.block<3,1>(0,ix+2) << x.block<3,1>(0,1) + Eigen::Vector3d(ho2otn, 0.0, -ho2);
    x.block<3,1>(0,ix+3) << x.block<3,1>(0,2) + Eigen::Vector3d(lmini, this->_width, 0.0);
    x.block<3,1>(0,ix+4) << x.block<3,1>(0,3) + Eigen::Vector3d(2.0 * ho2otn, 0.0, 2.0 * ho2);
    x.block<3,1>(0,ix+5) << x.block<3,1>(0,4) + Eigen::Vector3d(lmini, -this->_width, 0.0);
    x.block<3,1>(0,ix+6) << x.block<3,1>(0,5) + Eigen::Vector3d(ho2otn, 0.0, -ho2);
    x.block<3,1>(0,ix+7) << x.block<3,1>(0,6) + Eigen::Vector3d(pLength/2.0, 0.0, 0.0);
  //x.block<3,1>(0,ix+7) << Eigen::Vector3d(this->_totalLength / 2.0 , -this->_width, 0.0);

    x += pbcShift.replicate(1, 8);

    for (int i = ix; i < ix+8; ++i){
      Segment s(x.block<3,1>(0,i), x.block<3,1>(0,i+1), this->_burgers);
      subdivide(s);
    }
  }
  std::cout << "DEBUG: main nodes = " << std::endl;
  std::cout << x << std::endl;

}

void HelicalTurn::get_segments( std::vector<Segment> & storage) const
{
  storage.clear();
  storage = this->_segVector;
}

void HelicalTurn::write_to_file(FILE * pf) const
{
  fprintf ( pf, "#Helical turn %i segments, %1.4f total length",
            this->_segVector.size(), this->_totalLength);

  Eigen::Vector3d x;
  for (int i = 0; i < this->_segVector.size(); ++i){
    this->_segVector[i].get_head(x);
    fprintf (pf, "\n%1.5f %1.5f %1.5f", x(0), x(1), x(2));
  }

  this->_segVector.back().get_tail(x);
  fprintf (pf, "\n%1.5f %1.5f %1.5f", x(0), x(1), x(2));
}

double HelicalTurn::calculate_minimum_angle(double prismaticLength) const
{
  double tn = 2.0 * this->_height / (this->_totalLength - prismaticLength - 2.0 * this->_dx);
  return atan (tn);
}
