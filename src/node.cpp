#include "node.h"
#include "segment2.h"
#include "Eigen/Dense"
#include <vector>

Node::Node( Eigen::Vector3d& coords)
{
  _x << coords;
}

Node::~Node()
{
  //  if (_connectedSegments.size() > 0){
  //  printf ( "Error cannot delete Node %h, connected to %i segments",
  //           this, _connectedSegments.size() );
  //  return;
  // }
}

/*
Node* Node::connect(Segment2 *seg)
{
  _connectedSegments.push_back(seg);
  return this;
}

Node* Node::disconnect(Segment2 *seg)
{
  for (int i = 0; i < _connectedSegments.size(); ++i){
      if (seg == _connectedSegments[i]){
        _connectedSegments.erase(_connectedSegments.begin()+i);
        return this;
      }
  }
  printf ( "ERROR: trying to remove unconnected segment. Node = %p, Segment = %p\n", this, seg);
  return NULL;
}
*/

Eigen::Vector3d Node::get_coords() const
{
  return _x;
}

void Node::change_coords ( Eigen::Vector3d& newCoords )
{
  _x << newCoords;
}
