#include <iostream>
#include <vector>

auto sint(double& x, std::vector<double>& bex, std::vector<double>& rdbex,
  std::vector<double>& sex, std::vector<double>& betan, int b, double& alpha,
  double& wt, double& tbart, int nbx ){
  // Interpolates in scattering function, using SCT approximation to 
  // extrapolate outside the range in memory
  
  double sv, ex, sint, pi = 3.1415926535;
  if ( abs(x) > betan[b] ){
    if ( alpha <= 0.0 ){
      sv = 0.0;
    } else {
      ex = -(wt*alpha-abs(x))*(wt*alpha-abs(x))/(4*wt*alpha*tbart);
      if ( x > 0.0 ){ ex = ex - x; }
      sv = exp(ex/(4*pi*wt*alpha*tbart));
    }
    return sv;
  } 
  
  // interpolation
  int k1 = 1, k2 = b, k3 = nbx;
  // bisect for x
  int idone = 0;
  while ( idone == 0 ){
    if ( x == bex[k2] ){
      sv = sex[k2];
      return sv;
    }
    if ( x > bex[k2] ){
      k1 = k2;
      k2 = (k3-k2)/2 + k2;
      if ( k3-k1 <= 1 ){ idone = 1; }
    } 
    else {
      k3 = k2;
      k2 = (k2-k1)/2 + k1;
      if ( k3-k1 <= 1 ){ idone = 1; }
    }
  }
  std::cout << k1 << "    " << k2 << "     " << k3 << std::endl;
  
  std::cout << "Hello, world!" << std::endl;
  return 3.0;
}
