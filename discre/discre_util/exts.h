#include <iostream>
#include <vector>


auto exts( const std::vector<double>& sexpb, std::vector<double>& sex,
  const std::vector<double>& exb, const std::vector<double>& beta ){
  int k = beta.size();
  for ( auto b = 0; b < beta.size(); ++b ){
    sex[b] = sexpb[k-1];
    k -= 1;
  }
  if ( beta[0] <= 1.0e-9 ){
    sex[beta.size()-1] = sexpb[0];
    k = beta.size() + 1;
  }else{
    k = beta.size() + 2;
    sex[beta.size()] = sexpb[0];
  }
  for ( int b = 1; b < beta.size(); ++b ){
    sex[k-1] = sexpb[b]*exb[b]*exb[b];
    k += 1;
  }
  return 0;

}
