#include "segment.h"
#include "Eigen/Dense"
#include <math.h>

#define _USE_MATH_DEFINES

Segment::Segment(const Eigen::Vector3d& p1,
                 const Eigen::Vector3d& p2,
                 const Eigen::Vector3d& brg){
  this->_x1 << p1;
  this->_x2 << p2;
  this->_burgers << brg;

  this->_lineVector = (p2 - p1);
  this->_unitVector = this->_lineVector / this->_lineVector.norm();
}

double Segment::length() const{
  return this->_lineVector.norm();
}

void Segment::get_burgers_vector(Eigen::Vector3d & vec) const{
  vec << this->_burgers;
}

void Segment::get_unit_vector(Eigen::Vector3d &u) const {
  u << this->_unitVector;
}

void Segment::get_head(Eigen::Vector3d & u) const{
  u << this->_x1;
}

void Segment::get_tail(Eigen::Vector3d & u) const{
  u << this->_x2;
}

Segment::~Segment(){};
