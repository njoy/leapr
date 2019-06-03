#include <iostream>
#include <vector>

template <typename V, typename F>
auto addDeltaFuncs( const F twt, const F dwf, const V& bes, const V& betan, 
  const V& wts, V& sexpb, const int n ) {

  // Add the delta functions to the sexpb vec so it can be added to the 
  // scattering law 
  // delta(0.0) is saved for the incoherent elastic scattering
  
  if ( twt > 0.0  or dwf < 1.0e-10 ){ return; }

  unsigned int j; 
  int jj, m = 0;
  bool done = false;
  double add, db, be;

  while ( m < n ){
    m += 1;
    if ( bes[m-1] < 0.0 ){
      be = -bes[m-1];
      if ( be < betan[betan.size()-2] ){
        db = 1000;
        j = 0;
        while ( j < betan.size()-1 and not done ){
          j += 1;
          jj = j;
          if ( std::abs(be-betan[j-1] ) > db ){
            done = true;
          } else {
            db = std::abs(be-betan[j-1]);
          }
        }

        add = j > 2 ? 2 * dwf * wts[m-1]/(betan[jj-1]-betan[jj-3]) :
                          dwf * wts[m-1]/betan[jj-1];

        if ( add >= 1.0e-20 ){ sexpb[jj-2] += add; }
        if ( done ){ return; }
      }
    }
  }
}
