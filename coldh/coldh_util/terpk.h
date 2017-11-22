#include <iostream>
#include <vector>

auto terpk( std::vector<double>& ska, std::vector<double>& nka, double& delta,
  double& be ){
  /* Interpolates in a table of ska(kappa) for a required kappa
   */
  int i = be / delta;
  if ( i < nka - 1 ){
    double bt = i * delta;
    double btp = bt + delta;
    i += 1;
    return ska[i-1] + (be-bt)*(ska[i]-ska[i-1])/(delta);
  } else { 
    return 1.0;
  }
}

