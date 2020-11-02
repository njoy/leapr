#include <range/v3/all.hpp>

double factorial( int n ){
  return ( n <= 1 ) ? 1.0 : n*factorial(n-1);
}

auto getClebschGordon( int jj, int ll, int nn ){
  using std::pow;
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

  return pow(-1.0,(jj+ll-nn)/2)*sqrt((2.0*nn+1.0)/(jj+ll+nn+1))*c1*c2*c3*c4;
}



template <typename Float>
Float sjbes( int n, const Float& y ){
  /* Calculates spherical bessel functions for cold hydrogen or deuterium 
   * calculation. These will be used in Eq. 567 - 568
   * The spherical bessel function is of nth order, evaluated at position y.
   * In the case of cold hydrogen and deuterium calculations, y is defined by
   */
  int k, l, nm;
  Float w, sj, t1, t2, t3;
  // check for large or negative arguments
  if ( n >= 3.0e4 or y > 3.0e4 or y < 0.0 or n < 0.0 ){
    return 0.0;
  }
  // compute normal values
  if (y <= 7.0e-4) {
    if      (n  == 0) { return 1; } 
    else if (n  > 10) { return 0; } 
    else {
      t3 = 1;
      for ( int i = 3; i < 2*n + 3; i = i + 2 ){
        t3 = t3 * i;
      }
      return std::pow(y,n)/t3;
    }
  }
  else {
    w = y < 0.2 ? 1 - y*y * ( 1 - y*y/20 ) / 6 : 
                  sin(y) / y;
    if (n == 0) {
      return w;
    } 
    else {
      if      (y >= 100.0) { l = int( y/50 + 18 ); } 
      else if (y >= 10.0 ) { l = int( y/10 + 10 ); } 
      else if (y >  1.0  ) { l = int( y/2  + 5  ); } 
      else                 { l = 5; }

      nm = y > n ? y + l : n + l;
      t3 = 0;
      t2 = 2.0e-38;
      for ( auto i = 0; i < nm; ++i ){
        k = nm - i - 1;
        t1 = (2*k + 3) * t2 / y - t3;

        if (n == k){ sj = t1;}
        if (std::abs(t1) >= 1.0e25) {
          t1 = t1 / 1.0e25;
          t2 = t2 / 1.0e25;
          sj = sj / 1.0e25;
        }
        t3 = t2;
        t2 = t1;
      }
      return w * sj / t1;
    }
  }

}




template <typename Float>
auto sumh(int j, int jp, Float y){
  /* Does sum over Bessel functions and Clebsch-Gordon coefficients
   * for cold hydrogen or deuterium calculation, which follows Eq. 567 - 568
   */ 
  using std::pow;
  if      (j  == 0) { return pow(sjbes(jp,y) * getClebschGordon(j,jp,jp),2); } 
  else if (jp == 0) { return pow(sjbes(j ,y) * getClebschGordon(j,0, j ),2); }

  return ranges::accumulate(
           ranges::view::iota(std::abs(j-jp)+1,std::min(j+jp+2,std::abs(j-jp)+11)) 
         | ranges::view::transform( [y,j,jp](auto n){ 
             return pow(sjbes(n-1,y) * getClebschGordon(j,jp,n-1),2); }),
         0.0);

}
