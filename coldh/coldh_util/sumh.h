#include <iostream>
#include <vector>


auto sumh(int j, int jp, double y){
  /* Does sum over Bessel functions and Clebsch-Gordon coefficients
   * for cold hydrogen or deuterium calculation.
   */ 
  int imp, ipk1, mpk, ipk, n, n1;
  double sum1, sum2;
  if (j == 0) {
    sum2 = (sjbes(jp,y)*cn(j,jp,jp))*(sjbes(jp,y)*cn(j,jp,jp));
  } 
  else if (jp == 0) {
    sum2 = (sjbes(j,y)*cn(j,0,j))*(sjbes(j,y)*cn(j,0,j));
  }
  else {
    sum1 = 0;
    imk = iabs(j-jp) + 1;
    ipk1 = j + jp + 1;
    mpk = ipk1 - imk;
    if (mpk.le.9) {
      ipk = ipk1;
    } else {
      ipk = imk + 9;
    }
    for ( auto n = imk - 1; n < ipk; ++n ){
      n1 = n - 1;
      sum1 = sum1 + (sjbes(n1,y)*cn(j,jp,n1))*(sjbes(n1,y)*cn(j,jp,n1));
    }
    sum2 = sum1;
   }
   sumh = sum2;
   return;

}
