#include <iostream>
#include <vector>


void exts( const std::vector<double>& sexpb, std::vector<double>& sex,
  const std::vector<double>& exb, const std::vector<double>& beta ){

  // Here we reverse the vector sexpb and put it in the beginning
  // of sex.
  std::reverse_copy(std::begin(sexpb), std::end(sexpb), std::begin(sex) );

  int k = beta[0] <= 1.0e-9 ? beta.size() + 1 : beta.size() + 2;
  sex[k-2] = sexpb[0];

  for ( int b = 1; b < beta.size(); ++b ){
    sex[k-1] = sexpb[b]*exb[b]*exb[b];
    k += 1;
  }
}
