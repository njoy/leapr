#include <vector>
#include <range/v3/all.hpp>

double helper(int k, bool isEven ){
  if ( k <= 0 ){ return 1.0; }
  double s = ranges::accumulate(ranges::view::iota(1,k+1) 
           | ranges::view::transform([](auto x){ return log(x); }),0.0);
  std::cout << exp(s) << std::endl;
  double fact = (s > 0.0) ? exp(s) : 1.0;
  return isEven ? sqrt(fact) : fact;
}

auto cn( int jj, int ll, int nn ){
  /* Calculates Clebsch-Gordon coefficients for cold hydrogen or 
   * deuterium calculation. 
   */
  double c1, c2, c3, c4;

  std::cout << jj+ll+nn << "     " << helper( jj + ll + nn, true) << std::endl;
  // if sum of three inputs is even, continue. Else return 0.0
  if ( (jj + ll + nn )%2 == 0 ){

    c1 = helper(( jj + ll + nn ) * 0.5, false) / helper( jj + ll + nn, true);
    c2 = helper( jj + ll - nn, true) / helper(( jj + ll - nn ) * 0.5, false);
    c3 = helper( jj - ll + nn, true) / helper(( jj - ll + nn ) * 0.5, false);
    c4 = helper(-jj + ll + nn, true) / helper((-jj + ll + nn ) * 0.5, false);

    return std::pow(-1.0,(jj+ll-nn)/2) * sqrt((2.0*nn+1.0)/(jj+ll+nn+1)) * 
           c1 * c2 * c3 * c4;
  }
  else { return 0.0; }
}


