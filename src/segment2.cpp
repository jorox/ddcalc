#include "segment2.h"
#include "Eigen/Dense"
#include <math.h>

#define _USE_MATH_DEFINES

Segment2::Segment2(Node* p1,
                   Node* p2,
                   const Eigen::Vector3d& brg,
                   const Eigen::Matrix<int, 4, Eigen::Dynamic>& pln)
{
  //p1->connect(this);
  //p2->connect(this);
  this->_n1 = p1;
  this->_n2 = p2;
  this->_burgers << brg;

  _igplane = pln;

  this->refresh_data();
}

Segment2::Segment2(const Segment2& other)
{
  //other.head_node()->connect(this);
  //other.tail_node()->connect(this);

  this->_n1 = other.head_node();
  this->_n2 = other.tail_node();
  other.get_burgers_vector(this->_burgers);
  other.get_glide_plane(this->_igplane);

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

void Segment2::get_burgers_vector(Eigen::Vector3d & vec) const{
  vec << this->_burgers;
}

void Segment2::get_glide_plane(Eigen::Matrix<int, 4, Eigen::Dynamic> &plns) const
{
  plns = this->_igplane;
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
  //this->_n1->disconnect(this);
  //this->_n2->disconnect(this);
}
