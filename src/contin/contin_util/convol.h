#include <iostream>
template <typename F, typename A>
auto convol( const A& t1, const A& t2, const F& delta, A betaGrid,
  const int& nl=0, const int& n1=0, const int& nn=0 ){

  A t3(t2.size(),0.0);

  int i1, i2, len_t1 = int(t1.size()), len_t2 = int(t2.size());
  F f1, f2;
  len_t2 = nn;
  len_t1 = n1;

  for ( int i = 0; i < len_t2; ++i ){    // i iterates through t3
    for ( int j = 0; j < len_t1; ++j ){  // j iterates through t1, t2
      i1 = i + j;
      i2 = i - j;
      if ( t1[j] > 0 ){
        // Convolution will only be significant if t1[j] is not zero
    
        f1 = ( i1 - 1 > nl ) ? 0 : t2[i1]*exp(-betaGrid[j]);

        if      ( i2 >= 0 and  i2-1 < nl ){ f2 = t2[ i2];                   }
        else if ( i2 <  0 and -i2-3 < nl ){ f2 = t2[-i2] * exp( i2*delta ); }
        //else if ( i2 <  0 and -i2-3 < nl ){ f2 = t2[-i2] * exp( -betaGrid[i2] ); }
        else                              { f2 = 0;                         }

        // If at one of the endpoints, only give half contribution
        t3[i] += ( j == 0 or j == len_t1 - 1 ) ? t1[j] * ( f1 + f2 ) * 0.5 :
                                                 t1[j] * ( f1 + f2 );
      } // if
    } // for

    t3[i] = ( t3[i] * delta < 1e-30 ) ? 0 : t3[i] * delta;

  } // for
  return t3;

}

