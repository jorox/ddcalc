#include "node.h"
#include "segment2.h"
#include "faultedhelicalturn.h"
#include "Eigen/Dense"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

double lerp(double a, double b, double f)
{
  return a + f * (b - a);
}

int main( int argc, char** argv){

  double mu = 1.0;
  double nu = 0.2421;
  double cw = 1.4;

  //Cai ecalc(mu, nu, cw);
  double lmin = std::stod (argv[1]);
  double ltot = std::stod (argv[2]);
  double h = std::stod (argv[3]);
  double w = std::stod (argv[4]);
  double ang = std::stod (argv[5]) * M_PI / 180.;
  double plength = std::stod (argv[6]);

  FaultedHelicalTurn hturn(lmin, ltot, h, w, 1.0);

  FILE * out;
  out = fopen ( "faultedhturn.xml", "w");
  hturn.generate (plength, ang, 0);
  hturn.write_to_xml( out );



}
