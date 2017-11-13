#include <iostream>
#include <vector>


auto negativeTerms( int i, int& n, const double& normalizedEnergy, 
  const std::vector<double>& bminus, const int maxdd,
  std::vector<double>& wts, const std::vector<double>& wtn, 
  std::vector<double>& bes, const std::vector<double>& ben, int nn){

  int k = 0, idone = 0;
  double besn, wtsn;

  while ( k < 50 and idone == 0 ){
    k += 1;
    if ( bminus[k-1] <= 0 ){
      idone = 1;
    } else {
      for ( auto m = 0; m < nn; ++m ){
        besn = ben[m] - k * normalizedEnergy;
        wtsn = wtn[m] * bminus[k-1];
        if ( wtsn >= 1e-8 and n < maxdd ){
          n += 1;
          bes[n] = besn;
          wts[n] = wtsn;
        }
      }
    }
  }
  //for( auto entry : bes ){ std::cout << entry << std::endl; }
  //std::cout << wtsn << std::endl;
  return;
}
