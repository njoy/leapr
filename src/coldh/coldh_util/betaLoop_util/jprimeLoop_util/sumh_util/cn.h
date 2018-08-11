#include <vector>
#include <range/v3/all.hpp>

double factorial( int n ){
  if ( n <= 1 ){ return 1.0; }
  return n * factorial(n-1);
}

auto cn( int jj, int ll, int nn ){
  /* Calculates Clebsch-Gordon coefficients for cold hydrogen or 
   * deuterium calculation. 
   */
  double c1, c2, c3, c4;

  // if sum of three inputs is even, continue. Else return 0.0
  if ( (jj + ll + nn )%2 != 0 ){ return 0.0; }

  c1 =    factorial((jj+ll+nn)*0.5) / sqrt(factorial( jj+ll+nn));
  c2 = sqrt(factorial(( jj+ll-nn))) / factorial(( jj+ll-nn)*0.5);
  c3 = sqrt(factorial(( jj-ll+nn))) / factorial(( jj-ll+nn)*0.5);
  c4 = sqrt(factorial((-jj+ll+nn))) / factorial((-jj+ll+nn)*0.5);

  return std::pow(-1.0,(jj+ll-nn)/2) * sqrt((2.0*nn+1.0)/(jj+ll+nn+1)) * 
         c1 * c2 * c3 * c4;
}


