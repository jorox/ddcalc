#include "cai.h"
#include "helicalturn.h"
#include "segment.h"
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

  //double mu, lambda, nu, cw, gamma, _ALAT_, _MPANM2J_;

  //_ALAT_ = 0.323; //nm
  //_MPANM2J_ = 1e-21; // MJ/m^3 * nm^3 = 10^6 * 10^-27

  //mu = 36000; //MPa
  //lambda = 33800; //MPa
  //nu = lambda / 2. / (lambda + mu);
  //cw = 0.20; //nm
  //gamma = 0.135; // J/m^2

  //double N = 100;
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

  double mu = 1.0;
  double nu = 0.2421;
  double cw = 1.4;

  Cai ecalc(mu, nu, cw);
  double lmin = std::stod (argv[1]);
  double ltot = std::stod (argv[2]);
  double h = std::stod (argv[3]);
  double w = std::stod (argv[4]);
  double ang = std::stod (argv[5]) * M_PI / 180.;
  double plength = std::stod (argv[6]);

  HelicalTurn hturn(lmin, ltot, h, w, 1.0);

  std::vector<Segment> segs;
  double energy[2];

  int N = 100;
  double minang = hturn.calculate_minimum_angle (plength);
  double dlprism = (plength) / N;

  FILE * eturn;
  FILE * out;
  eturn = fopen ( "eturn.dat", "w");
  std::string fname;
  double hturnLength;

  double lprism;
  for (int i = 0; i < N; ++i){
    lprism = i * dlprism;
    hturn.regenerate(lprism, ang, 0);
    hturn.get_segments(segs);
    ecalc.calculate_total_energy(segs, energy);
    hturnLength = hturn.get_total_length();

    printf ( "lprism = %f, self-energy = %e, inter-energy = %e, length = %f\n",
             lprism, energy[0],  energy[1], hturnLength);

    fprintf ( eturn, "%f %e %e %e\n", lprism,
              energy[0] / hturnLength,  energy[1] / hturnLength,
              (energy[0] + energy[1]) / hturnLength );

    fname = std::string("hturn") + std::to_string(i) + std::string(".dat");
    out  = fopen(fname.c_str(), "w");

    hturn.write_to_file(out);
    fclose(out);
  }
  fclose( eturn );


}
