#include <iostream>
#include <vector>

auto terpk( std::vector<double>& ska, double& delta, double& x ){
  /* Interpolates in a table of ska(kappa) for a required kappa
   */
  int i = x / delta;
  return i < ska.size() - 1 ? ska[i] + (x-i*delta) * (ska[i+1]-ska[i]) / delta :
                              1.0;
}

