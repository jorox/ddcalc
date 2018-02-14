#include "cai.h"
#include "segment.h"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <math.h>
#define _USE_MATH_DEFINES

Cai::Cai( double mu, double nu, double cw){
  this->_mu = mu;
  this->_nu = nu;
  this->_cw = cw;
}

double Cai::special_func(const Eigen::Vector3d& R,
                         const Eigen::Vector3d& t1,
                         const Eigen::Vector3d& t2,
                         const Eigen::Vector3d& b1,
                         const Eigen::Vector3d& b2) const
{
  double b1b2 = b1.dot(b2);
  double t1t2 = t1.dot(t2);
  double b1t1 = b1.dot(t1);
  double b2t2 = b2.dot(t2);
  double b1t2 = b1.dot(t2);
  double b2t1 = b2.dot(t1);

  Eigen::Vector3d u = t1.cross(t2);
  Eigen::Vector3d v = u.cross(t1);
  Eigen::Vector3d vp = t2.cross(u);
  double Ra = sqrt(R.squaredNorm() + _cw * _cw);
  double W0 = _mu / 4. / M_PI / (1 - _nu) / u.dot(u);

  double uu = u.dot(u);
  double Ru = R.dot(u);
  double Rv = R.dot(v);
  double Rvp = R.dot(vp);
  double Rt1 = R.dot(t1);
  double Rt2 = R.dot(t2);

  double A1 = ( (1-_nu) * b1t1 * b2t2 + 2. * _nu * b2t1 * b1t2 );
  double A2 = ( b1b2 + b1t1 * b2t2 ) * t1t2;
  double A2p= ( b1b2 + b1t2 * b2t1 ) * t1t2;
  double A3p= ( (b1.dot(u) * b2.dot(vp) + b2.dot(u) * b1.dot(vp)) *
                t1t2 / uu );
  double A3 = ( (b1.dot(u) * b2.dot(v) + b2.dot(u) * b1.dot(v)) *
                t1t2 / uu );
  double A4 = ( b1t1 * b2.dot(v) + b1t2 * b2.dot(vp) ) * t1t2;
  double A5 =  b1.cross(u).dot( b2.cross(u) ) * 2.0 * t1t2 / uu;

  double res = 0;
  res +=  (A1 - A2p) * Rvp * log(Ra + Rt2);
  res += A3p * Ru * log(Ra + Rt2);
  res += (A1 - A2) * Rv * log(Ra + Rt1);
  res += A3 * Ru * log(Ra + Rt1);
  res += A4 * Ra;
  res += ( (A1 - A5) * (2.0 * Ru*Ru + uu * _cw*_cw) /
           sqrt(uu * _cw*_cw + Ru*Ru) *
           atan( ((1. + t1t2) * Ra + R.dot(t1+t2))/
                 sqrt(uu * _cw*_cw + Ru*Ru) )
           );

  return res * W0;
}

double Cai::special_func_parallel( const Eigen::Vector3d& R,
                                   const Eigen::Vector3d& t,
                                   const Eigen::Vector3d& b1,
                                   const Eigen::Vector3d& b2) const
{

  double b1DotT = b1.dot(t);
  double b1DotR = b1.dot(R);
  double b2DotT = b2.dot(t);
  double b2DotR = b2.dot(R);
  double b1DotB2 = b1.dot(b2);
  double tDotR = R.dot(t);
  double Ra = sqrt( R.squaredNorm() + _cw * _cw);
  double res = 0.;

  res = ( b1DotT * b2DotR + b1DotR * b2DotT -
          ( (2-_nu) * b1DotT * b2DotT + b1DotB2 ) * tDotR );
  res *= log(Ra + tDotR);
  res += ( (1-_nu) * b1DotT * b2DotT + b1DotB2 ) * Ra;
  res -= ( ( b1DotR - tDotR * b1DotT ) * ( b2DotR - tDotR * b2DotT ) / (Ra * Ra - tDotR * tDotR) ) * Ra;
  res += _cw * _cw *( ( (1+_nu) * b1DotT * b2DotT - 2.0 * b1DotB2 ) / ( 2.0 * (Ra * Ra - tDotR * tDotR) ) ) * Ra;
  res *= _mu / 4.0 / M_PI / (1 - _nu);
  return res;
}

double Cai::self_energy( const Segment& neo)
{
  double L = neo.length();
  Eigen::Vector3d burg, t;
  neo.get_burgers_vector(burg);
  neo.get_unit_vector(t);

  double bDotB = burg.dot(burg);
  double bDotT = burg.dot(t);
  double La = sqrt(L * L + _cw * _cw);

  double wSelf  = (bDotB - _nu * bDotT * bDotT) * L * log( (La + L)/_cw );
  wSelf -= (3-_nu) / 2.0 * bDotT * bDotT * (La - _cw);

  wSelf *= _mu / 4.0 / M_PI / (1-_nu);
  return wSelf;
}

double Cai::interaction_energy(const Segment& neo,
                               const Segment& other)
{
  Eigen::Vector3d t1, t2, b1, b2, p1, p2, p3, p4;
  Eigen::Vector3d R24, R14, R13, R23;
  double L24, L14, L13, L23;

  neo.get_unit_vector(t1);
  other.get_unit_vector(t2);
  neo.get_burgers_vector(b1);
  other.get_burgers_vector(b2);
  neo.get_head(p1);
  neo.get_tail(p2);
  other.get_head(p3);
  other.get_tail(p4);

  R24 = p4 - p2;
  L24 = R24.norm();
  R13 = p3 - p1;
  L13 = R13.norm();
  R14 = p4 - p1;
  L14 = R14.norm();
  R23 = p3 - p2;
  L23 = R23.norm();

  double wInt=0;
  if ( fabs(t1.dot(t2)) > 0.995 ){
    wInt += special_func_parallel(R24, t1, b1, b2);
    wInt += special_func_parallel(R13, t1, b1, b2);
    wInt -= special_func_parallel(R14, t1, b1, b2);
    wInt -= special_func_parallel(R23, t1, b1, b2);
  }
  else{
    wInt += special_func(R24, t1, t2, b1, b2);
    wInt += special_func(R13, t1, t2, b1, b2);
    wInt -= special_func(R14, t1, t2, b1, b2);
    wInt -= special_func(R23, t1, t2, b1, b2);
  }

  return wInt;
}

void Cai::calculate_total_energy ( const std::vector<Segment>& segs,
                                   double* energy)
{
  energy[0] = 0.0;
  energy[1] = 0.0;
  double wint=0.0;

  for (int i = 0; i < segs.size(); ++i){
    energy[0] += self_energy(segs[i]);
    wint = 0.0;
    for (int j = i+1; j < segs.size(); j++){
      wint += interaction_energy (segs[i], segs[j]);
    }
    energy[1] += wint;
  }

};
