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

  //this->generate(this->_totalLength / 2.0, M_PI / 2.0, 0);
  //FILE * fout;
  //fout = fopen("construct.dat","w");
  //  this->write_to_file(fout);
  //fclose(fout);
}

void HelicalTurn::subdivide(const Segment &s1)
{
  double dx = this->_dx;
  int n = (int) (s1.length() / dx);
  double rem = (s1.length() - (double) (n * dx)); // remainder can be zero
  Eigen::Vector3d start;
  Eigen::Vector3d end;
  Eigen::Vector3d t;
  s1.get_head(start);
  s1.get_tail(end);
  s1.get_unit_vector(t);

  // psuh first segment
  Segment tmp(start, start + t * dx, this->_burgers);
  this->_segVector.push_back(tmp);
  const Segment * p0 = &(this->_segVector.back());
  //printf( "DEBUG: start = %f %f %f\n", start(0), start (1), start(2));

  // push middle segments
  for (int i = 1; i < n; ++i){
    this->_segVector.back().get_tail(start);
    Segment tmp(start, start + t * dx, this->_burgers);
    this->_segVector.push_back(tmp);
    //printf( "DEBUG: start = %f %f %f\n", start(0), start (1), start(2));
  }

  // push last segment which will probably have length < dx
  this->_segVector.back().get_tail(start);
  //printf( "DEBUG: start = %f %f %f\n", start(0), start (1), start(2));
  if (rem > 0.001){
    Segment tmp(start, start + rem * t, this->_burgers);
    this->_segVector.push_back(tmp);
  }

  p0->get_head(start);
  this->_segVector.back().get_tail(end);

  printf ("   generated %i segments, (%f, %f, %f) --> (%f, %f, %f)\n",
          n, start(0), start(1), start(2), end(0), end(1), end(2));
}

void HelicalTurn::generate ( double pLength,
                             double inclination,
                             int image = 0)
{
  if (pLength < this->_dx){ pLength = 0.; }
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
  Eigen::Matrix<double, 3, Eigen::Dynamic> x(3, nMainNodes);
  Eigen::Matrix<double, 3, 1> pbcShift;
  int ix = 0;

  for (int ipbc = -image; ipbc < image + 1; ++ipbc){
    ix = ipbc + image;
    pbcShift << double(ipbc) * this->_totalLength, 0.0, 0.0;

    x.block<3,1>(0,ix) << -this->_totalLength / 2.0 , -this->_width, 0.0;
    x.block<3,1>(0,ix+1) << x.block<3,1>(0,ix) + Eigen::Vector3d(pLength/2.0, 0.0, 0.0);
    x.block<3,1>(0,ix+2) << x.block<3,1>(0,ix+1) + Eigen::Vector3d(ho2otn, 0.0, -ho2);
    x.block<3,1>(0,ix+3) << x.block<3,1>(0,ix+2) + Eigen::Vector3d(lmini, this->_width, 0.0);
    x.block<3,1>(0,ix+4) << x.block<3,1>(0,ix+3) + Eigen::Vector3d(2.0 * ho2otn, 0.0, 2.0 * ho2);
    x.block<3,1>(0,ix+5) << x.block<3,1>(0,ix+4) + Eigen::Vector3d(lmini, -this->_width, 0.0);
    x.block<3,1>(0,ix+6) << x.block<3,1>(0,ix+5) + Eigen::Vector3d(ho2otn, 0.0, -ho2);
    x.block<3,1>(0,ix+7) << x.block<3,1>(0,ix+6) + Eigen::Vector3d(pLength/2.0-this->_dx, 0.0, 0.0);

    x.block<3, 8>(0, ix) += pbcShift.replicate(1, 8);
  }
  std::cout << "DEBUG: main nodes = " << std::endl;
  std::cout << x << std::endl;
  std::cout << "DEBUG: subdividing ... " << std::endl;

  for (int i = 0; i < x.cols()-1 ; ++i){
    Segment s (x.block<3,1>(0,i), x.block<3,1>(0,i+1), this->_burgers);
    if (s.length() >= this->_dx){
      subdivide(s);
    }
  }
  std::cout << "DEBUG: subdivision complete " << std::endl;
}
/*
void HelicalTurn::get_segments( std::vector<Segment> & storage) const
{
  storage.clear();
  for (int i = 0; i < this->_segVector.size(); ++i){
    storage.push_back(Segment(this->_segVector[i]));
  }
}
*/
void HelicalTurn::write_to_file(FILE * pf) const
{
  fprintf ( pf, "#Helical turn %i segments, %1.4f total length",
            this->_segVector.size(), this->_totalLength);

  Eigen::Vector3d x;
  Eigen::Vector3d b;
  for (int i = 0; i < this->_segVector.size(); ++i){
    this->_segVector[i].get_head(x);
    this->_segVector[i].get_burgers_vector(b);
    fprintf (pf, "\n%1.5f %1.5f %1.5f %1.5f", x(0), x(1), x(2), b.norm());
  }

  this->_segVector.back().get_tail(x);
  this->_segVector.back().get_burgers_vector(b);
  fprintf (pf, "\n%1.5f %1.5f %1.5f %1.5f", x(0), x(1), x(2), b.norm());
}

double HelicalTurn::calculate_minimum_angle(double prismaticLength) const
{
  double tn = 2.0 * this->_height / (this->_totalLength - prismaticLength - 2.0 * this->_dx);
  return atan (tn);
}

double HelicalTurn::get_total_length() const
{
  double res = 0;
  for (int i = 0; i < _segVector.size(); ++i){
    res += _segVector[i].length();
  }
  return res;
}
