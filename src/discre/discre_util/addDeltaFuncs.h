#include <iostream>
#include <vector>

template <typename Float, typename Range>
auto addDeltaFuncs( const Float twt, const Float& dwf, const Range& bes, 
  const Range& betan, const Range& wts, Range& sexpb, const int n ) {
  // Add the delta functions to the scattering law 
  // delta(0.0) is saved for the incoherent elastic scattering
  
  if ( twt > 0.0  or dwf < 1.0e-10 ){ return; }

  unsigned int j; 
  int jj, m = 0;
  bool done = false;
  Float add, db;

  while ( m < n ){
    m += 1;
    if ( bes[m-1] < 0.0 ){
      if ( -bes[m-1]< betan[betan.size()-2] ){
        db = 1000;
        j = 0;
        while ( j < betan.size()-1 and not done ){
          j += 1;  //jj = j;
          if ( std::abs(-bes[m-1]-betan[j-1] ) > db ){
            done = true;
          } else {
            db = std::abs(-bes[m-1]-betan[j-1]);
          }
        }
        jj = j;

        add = j > 2 ? 2 * dwf * wts[m-1]/(betan[jj-1]-betan[jj-3]) :
                          dwf * wts[m-1]/betan[jj-1];

        if ( add >= 1.0e-20 ){ sexpb[jj-2] += add; }
        if ( done ){ return; }
      }
    }
  }
}
