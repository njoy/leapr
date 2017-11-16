#include <iostream>
#include <vector>

auto sint(const double& x, const std::vector<double>& bex, 
  const std::vector<double>& rdbex, const std::vector<double>& sex,
  const std::vector<double>& betan, int b, const double& alpha,
  const double& wt, const double& tbart, int nbx ){
  // Interpolates in scattering function, using SCT approximation to 
  // extrapolate outside the range in memory
  
  std::cout << "\n" << std::endl;
  double sv, ex, sint, pi = 3.1415926;
  if ( abs(x) > betan[b] ){
    if ( alpha <= 0.0 ){                     // OPTION A
      std::cout <<"A"<<std::endl;    
      sv = 0.0;
    } else {                                 // OPTION B1
      std::cout <<"B"<<std::endl;
      ex = -(wt*alpha-abs(x))*(wt*alpha-abs(x))/(4*wt*alpha*tbart);
      if ( x > 0.0 ){ 
        std::cout <<"B2"<<std::endl;        // OPTION B2
        ex = ex - x; 
      }
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
      std::cout <<"C"<<std::endl;
      return sv;
    }
    if ( x > bex[k2] ){               
      k1 = k2;
      k2 = (k3-k2)/2 + k2;
      if ( k3-k1 <= 1 ){                   // OPTION D
        idone = 1; 
        std::cout <<"D"<<std::endl;
      }
    } 
    else {
      k3 = k2;
      k2 = (k2-k1)/2 + k1;                 // OPTION E
      if ( k3-k1 <= 1 ){ 
        idone = 1;
        std::cout <<"E"<<std::endl;
      }
    }
  }
  std::cout << k1 << "    " << k2 << "     " << k3 << std::endl;
  double ss1, ss3;
  if ( sex[k1] <= 0.0 ){
    ss1 = -225.0;                         // OPTION F
    std::cout <<"F"<<std::endl;
  } else {
    ss1 = log( sex[k1] );                 // OPTION G
    std::cout << sex[k1-1] << std::endl;
    std::cout <<"G"<<std::endl;
  }
  if ( sex[k3] <= 0.0 ){
    ss3 = -225.0;                         // OPTION H
    std::cout <<"H"<<std::endl;
  } else {
    ss3 = log( sex[k3] );                 // OPTION I
    std::cout <<"I"<<std::endl;
  }
  std::cout << "ss1:   " << ss1 << std::endl;
  std::cout << "ss3:   " << ss3 << std::endl;
  
  ex = ( (bex[k3]-x)*ss1+(x-bex[k1])*ss3 ) * rdbex[k1];
  std::cout << "ex:    " << ex << std::endl;
  std::cout << rdbex[k1] << std::endl;
  sv = 0;
  if ( ex > -225.0 ){ 
    sv = exp(ex);                        // OPTION J
    std::cout <<"J"<<std::endl;
  }
  std::cout << "sv:    " << sv << std::endl;
  return sv;

  


  
  std::cout << "Hello, world!" << std::endl;
  return 3.0;
}
