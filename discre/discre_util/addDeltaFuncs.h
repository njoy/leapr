#include <iostream>
#include <vector>

auto addDeltaFuncs( double twt, double dwf, std::vector<double>& bes,
    std::vector<double>& betan, std::vector<double>& wts, 
    std::vector<double>& sexpb, int n) {

  // Add the delta functions to the scattering law 
  // delta(0.0) is saved fro the incoherent elastic scattering
  int j, jj;
  double add;
    std::cout << twt << std::endl;
  if ( twt <= 0.0 ){
    std::cout << "HERE" << std::endl;
    int m = 0, idone = 0;
    while ( m < n and idone == 0 ){
      std::cout << m << std::endl;
      m += 1;
      if ( dwf < 1.0e-10 ){
        idone = 1;
      } else {
        if ( bes[m-1] < 0.0 ){
          double be = -bes[m-1];
          if ( be < betan[betan.size()-2] ){
            double db = 1000;
            idone = 0;
            j = 0;
            while ( j < betan.size()-1 and idone == 0 ){
              j += 1;
              int jj = j;
              if ( abs(be-betan[j-1] )>db ){
                idone = 1;
              } else {
                db = abs(be-betan[j-1]);
              }
            }
            if ( jj <=2 ){
              add = wts[m-1]/betan[jj-1];
            } else {
              add = 2 * wts[m-1]/(betan[jj-1]-betan[jj-2]);
            }
            add *= dwf;
            if ( add >= 1.0e-20 ){ sexpb[jj-2] += add; }
          }
        }
      }
    }
  }
      

  std::cout << "Hello, world!" << std::endl;
}
