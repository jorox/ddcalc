#include "faultedhelicalturn.h"
#include "node.h"
#include "segment2.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#define _USE_MATH_DEFINES

FaultedHelicalTurn::FaultedHelicalTurn( double lmin,
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
}

void FaultedHelicalTurn::generate ( double pLength,
                                    double inclination,
                                    int period = 0)
{
  double minInc = this->calculate_minimum_angle(pLength);
  printf( "prismatic segment = %f, min angle = %f\n", pLength, minInc * 180. / M_PI);
  if (inclination < minInc) {
    printf ( "ERROR: FaultedHelicalTurn::regenerate, minimum possible angle = %f degrees\n", minInc * 180. / M_PI);
    return;
  }

  //_segVector.clear();

  double tn = tan( inclination );
  double ho2 = this->_height / 2.0;
  double ltop = ( _totalLength - pLength - 2.0 * _height / tn ) / 2.0;
  double d = _height / 2.0 / tn;
  printf ("DEBUG: l = %f\n", ltop);


  // generate first prismatic segment
  int nMainNodes = 6 + 2 * 6 * period;
  Eigen::Matrix<double, 3, Eigen::Dynamic> x(3, nMainNodes);
  Eigen::Matrix<double, 3, 1> pbcShift;
  int ix = 0;

  for (int ipbc = -period; ipbc < period + 1; ++ipbc){
    ix = ipbc + period;
    pbcShift << double(ipbc) * this->_totalLength, 0.0, 0.0;

    x.block<3,1>(0,ix) << -2.0 * d - ltop , -_width, 0;
    x.block<3,1>(0,ix+1) << -d - ltop, -_width, ho2;
    x.block<3,1>(0,ix+2) << -d,  0., ho2;
    x.block<3,1>(0,ix+3) <<  d, 0., -ho2;
    x.block<3,1>(0,ix+4) <<  d + ltop, -_width, -ho2;
    x.block<3,1>(0,ix+5) <<  2.0 * d + ltop, -_width, 0.;

    x.block<3, 6>(0, ix) += pbcShift.replicate(1, 6);
  }

  Eigen::Vector3d p;
  for (int i = 0; i < x.cols(); ++i){
    p << x.block<3,1>(0,i);
    _nodeVector.push_back(Node(p));
  }

  Eigen::Matrix<int,4,Eigen::Dynamic> gprism(4,1);
  Eigen::Matrix<int,4,Eigen::Dynamic> gbasal(4,1);
  gprism << -1,1,0,0;
  gbasal << 0,0,0,1;

  //push first prismatic across periodic boundary
  _segVector.push_back(Segment2(&_nodeVector.back(), &_nodeVector.front(),
                                _burgers, gprism));
  for (int i = 0; i < x.cols()-1; ++i){
    if (i%2 == 0){
      _segVector.push_back(Segment2(&_nodeVector[i], &_nodeVector[i+1],
                                    _burgers, gbasal));
    }
    else{
      _segVector.push_back(Segment2(&_nodeVector[i], &_nodeVector[i+1],
                                    _burgers, gprism));
    }
  }
  printf ("... done generating %i segments\n", _segVector.size());
}

/**
void FaultedHelicalTurn::write_to_file(FILE * pf) const
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
**/

double FaultedHelicalTurn::calculate_minimum_angle(double prismaticLength) const
{
  double tn = 2.0 * this->_height / (this->_totalLength - prismaticLength - 2.0 * this->_dx);
  return atan (tn);
}

double FaultedHelicalTurn::get_total_length()
{
  double res = 0;
  for (int i = 0; i < _segVector.size(); ++i){
    res += _segVector[i].length();
  }
  return res;
}

int FaultedHelicalTurn::find_node_tag(Node *targetNode)
{
  for (int i = 0; i < _nodeVector.size(); ++i){
    if (targetNode == &_nodeVector[i]) return i;
  }
  return -1;
}

void FaultedHelicalTurn::write_to_xml(FILE *pf)
{
  int inode;
  Eigen::Vector3d xnode;
  fprintf (pf, "<?xml version=\"1.0\" ?>\n" );
  fprintf (pf, "<root>\n\t<nodes>");
  for (int i = 0; i < _nodeVector.size(); ++i){
    xnode = _nodeVector[i].get_coords();
    fprintf (pf, "\n\t\t<node tag=\"%i\" pinned=\"0\" x=\"%1.8f\" y=\"%1.8f\" z=\"%1.8f\" />",
             i, xnode(0), xnode(1), xnode(2));
  }
  fprintf (pf, "\n\t</nodes>");
  fprintf (pf, "\n\t<lines>");
  Eigen::Matrix<int,4,Eigen::Dynamic> gplanes;
  for (int i = 0; i < _segVector.size(); ++i){
    _segVector[i].get_glide_plane(gplanes);
    fprintf (pf, "\n\t\t<line tag=\"%i\" bh=\"2\" bk=\"2\" bi=\"-4\" bl=\"0\">", i);
    fprintf (pf, "\n\t\t\t<constraints>");
    for (int j = 0; j < gplanes.cols(); ++j){
      fprintf (pf, "\n\t\t\t\t<plane h=\"%i\" k=\"%i\" i=\"%i\" l=\"%i\" />",
               gplanes(0,j), gplanes(1,j), gplanes(2,j), gplanes(3,j));
    }
    fprintf (pf, "\n\t\t\t</constraints>");
    fprintf (pf, "\n\t\t\t<nodes>");
    fprintf (pf, "\n\t\t\t\t<node tag=\"%i\" />",
             find_node_tag(_segVector[i].head_node()));
    fprintf (pf, "\n\t\t\t\t<node tag=\"%i\" />",
             find_node_tag(_segVector[i].tail_node()));
    fprintf (pf, "\n\t\t\t</nodes>");
    fprintf (pf, "\n\t\t</line>");
  }
  fprintf (pf, "\n\t</lines>");
  fprintf (pf, "\n\t<entangled_nodes />");
  fprintf (pf, "\n</root>");
}
