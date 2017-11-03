#include <iostream>
#include <vector>
#include <cmath>

auto convol( const std::vector<double>& t1, const std::vector<double>& t2, 
             const double& delta ){

  std::vector<double> t3(t2.size(),0.0);

  int i1, i2;
  double f1, f2;

  for ( int i = 0; i < t2.size(); ++i ){    // i iterates through t3
    for ( int j = 0; j < t1.size(); ++j ){  // j iterates through t1, t2
      i1 = i + j;
      i2 = i - j;
      if ( t1[j] > 0 ){   // continue if the jth entry of t1 is positive
                          // If it's zero or negative, ignore since
                          // it won't contribute to to t3
        f1 = 0.0;
        if ( i1 - 1 <=  int(t1.size()) ){
          f1 = t2[ i1 ]*exp( -j*delta );
        }

        f2 = 0;
        if ( i2 >= 0 and   i2  <= int(t1.size()) ){ f2 = t2[ i2 ]; }
        if ( i2 <  0 and -i2-2 <= int(t1.size()) ){ f2 = t2[-i2 ] *
                                                          exp( i2 * delta ); };

        if ( 0 < j and j < t1.size() - 1 ){
          t3[i] += t1[j] * (f1+f2);
        } else { 
          t3[i] += t1[j] * (f1+f2) * 0.5;
        }

      } // if
    } // for

    t3[i] = ( t3[i] * delta < 1e-30 ) ? 0 : t3[i] * delta;

  } // for
  return t3;
}


