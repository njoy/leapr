#include <iostream>
#include <vector>

auto sint(const double& x, const std::vector<double>& bex, 
  const std::vector<double>& rdbex, const std::vector<double>& sex,
  const std::vector<double>& betan, int b, const double& alpha,
  const double& wt, const double& tbart, int nbx ){
  // Interpolates in scattering function, using SCT approximation to 
  // extrapolate outside the range in memory
  
  double sv, ex, sint, pi = 3.1415926;
  if ( abs(x) > betan[b] ){
    if ( alpha <= 0.0 ){                     // OPTION A
      sv = 0.0;
    } else {                                 // OPTION B
      ex = -(wt*alpha-abs(x))*(wt*alpha-abs(x))/(4*wt*alpha*tbart);
      if ( x > 0.0 ){ ex = ex - x; }
      sv = exp(ex)/(4*pi*wt*alpha*tbart);
    }
    return sv;
  } 
  
  // interpolation
  int k1 = 1, k2 = b, k3 = nbx;
  // bisect for x
  int idone = 0;
  while ( idone == 0 ){
    if ( x == bex[k2] ){                   // OPTION C
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
