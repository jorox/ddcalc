#include "cai.h"
#include "helicalturn.h"
#include "segment.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

double lerp(double a, double b, double f)
{
  return a + f * (b - a);
}

int main( int argc, char** argv){

  double mu, lambda, nu, cw, gamma, _ALAT_, _MPANM2J_;

  _ALAT_ = 0.323; //nm
  _MPANM2J_ = 1e-21; // MJ/m^3 * nm^3 = 10^6 * 10^-27

  mu = 36000; //MPa
  lambda = 33800; //MPa
  nu = lambda / 2. / (lambda + mu);
  cw = 0.20; //nm
  gamma = 0.135; // J/m^2

  double N = 100;
  /**
  double maxSep = 5.0; //nm
  double minSep = 0.2; //nm
  double step = (maxSep - minSep) / N;
  double sep = minSep;
  Segment s1 ( Eigen::Vector3d(-5.0, 0.0, 0.0),
               Eigen::Vector3d(5.0, 0.0, 0.0),
               Eigen::Vector3d(_ALAT_/2.0, 0.0, 0.0));
  Segment s2 ( Eigen::Vector3d(-5.0, sep, 0.0),
               Eigen::Vector3d(5.0, sep, 0.0),
               Eigen::Vector3d(_ALAT_/2.0, 0.0, 0.0) );

  std::ofstream fout;
  fout.open("result.txt");
  fout << "#distance(nm) energy(J/m)\n";

  for (int i = 0; i<N; ++i){
    Segment s2 ( Eigen::Vector3d(-5.0, sep, 0.0),
                 Eigen::Vector3d(5.0, sep, 0.0),
                 Eigen::Vector3d(_ALAT_/2.0, 0.0, 0.0) );

    double intW = cai::interaction_energy(mu, nu, cw, s1, s2);
    double slfW = cai::self_energy(mu, nu, cw, s1) + cai::self_energy(mu, nu, cw, s2);

    std::cout << "interaction energy (J) = ";
    std::cout << intW * _MPANM2J_;
    std::cout << std::endl;
    std::cout << "self-energy (J)= ";
    std::cout << slfW * _MPANM2J_;
    std::cout << std::endl;

    fout << sep << " " << intW * _MPANM2J_ + slfW * _MPANM2J_ + gamma * sep * s1.length() * 1e-18 << std::endl;

    sep += step;
  }
  fout.close();
  return 0;
  **/

  double lmin = std::stod (argv[1]);
  double ltot = std::stod (argv[2]);
  double h = std::stod (argv[3]);
  double w = std::stod (argv[4]);
  double ang = std::stod (argv[5]) * M_PI / 180.0;
  printf( "DEBUB ang = %f\n", ang);
  double plength = std::stod (argv[6]);
  HelicalTurn hturn(lmin, ltot, h, w, 3.232);
  hturn.regenerate(plength, ang);
  FILE * out;
  out  = fopen("hturn.dat", "w");

  hturn.write_to_file(out);


}
