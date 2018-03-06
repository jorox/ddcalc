#include "segment2.h"
#include "Eigen/Dense"
#include <math.h>

#define _USE_MATH_DEFINES

Segment2::Segment2(Node* p1,
                   Node* p2,
                   const Eigen::Vector3d& brg)
{
  p1->connect(this);
  p2->connect(this);
  this->_n1 = p1;
  this->_n2 = p2;
  this->_burgers << brg;

  this->refresh_data();
}

void Segment2::refresh_data()
{
  this->_x1 << _n1->get_coords();
  this->_x2 << _n2->get_coords();
  this->_lineVector = (this->_x2 - this->_x1);
  this->_unitVector = this->_lineVector / this->_lineVector.norm();
}

double Segment2::length(){
  this->refresh_data();
  return this->_lineVector.norm();
}

void Segment2::get_burgers_vector(Eigen::Vector3d & vec){
  this->refresh_data();
  vec << this->_burgers;
}

Node* Segment2::head_node() const
{
  return this->_n1;
}

Node* Segment2::tail_node() const
{
  return this->_n2;
}

Segment2::~Segment2()
{
  this->_n1->disconnect(this);
  this->_n2->disconnect(this);
}
