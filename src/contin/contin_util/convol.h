#include <range/v3/all.hpp>

template <typename Float, typename Range>
auto getConvolAtPoint( int i, const Float& delta, const Range& t1, 
  const Range& t2 ){
  Float sumVal = 0.0, toAdd, expVal;
  int len_t1 = t1.size();
  for ( int j = -len_t1+1; j < len_t1; ++j ){  // j iterates through t1, t2
    if ( i - j >= int(t2.size()) ){ return sumVal; }
    expVal = j < 0 ? exp(j*delta) :
        ( i < j ? exp((i-j)*delta) :
          1.0 );
    toAdd = t1[abs(j)]*t2[abs(i-j)]*expVal;
    sumVal += ( j == -len_t1+1 or j == len_t1-1 ) ? 0.5*toAdd : toAdd;
  } 
  return sumVal;
}


template <typename Float, typename Range>
auto convol( const Range& t1, const Range& t2, const Float& delta, const int nn){
  Range t3(nn,0.0);
  for ( int i = 0; i < nn; ++i ){   
    t3[i] = getConvolAtPoint(i,delta,t1,t2) * delta;
    t3[i] = ( t3[i] < 1e-30 ) ? 0 : t3[i];
  } 
  return t3;
}


/*
template <typename F, typename A>
auto getConvolAtPoint( int i, int nl, F delta, const A& t1, const A& t2, int len_t1 ){
    F sumVal = 0.0;
    for ( int j = 0; j < len_t1; ++j ){  // j iterates through t1, t2
      if ( t1[j] > 0 ){
        // Convolution will only be significant if t1[j] is not zero
        F f1, f2;
        //---------------------------------------------------------------------
        int i1 = i + j;
        f1 = ( i1 < nl + 2 ) ? t2[i1]*exp(-j*delta) : 0;
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
        int i2 = i - j;
        if      ( i2 >= 0 and     i2  < nl+1 ){ f2 = t2[ i2];                   }
        else if ( i2 <  0 and abs(i2) < nl+3 ){ f2 = t2[-i2] * exp( i2*delta ); }
        else                                  { f2 = 0;                         }
        //---------------------------------------------------------------------
        // If at one of the endpoints, only give half contribution
        sumVal += ( j == 0 or j == len_t1 - 1 ) ? t1[j] * ( f1 + f2 ) * 0.5 :
                                                  t1[j] * ( f1 + f2 );
      } // if
    } // for
    return sumVal;
}
*/


/*
template <typename F, typename A>
auto convolOriginal( const A& t1, const A& t2, const F& delta, const int& nl=0,
  const int& nn=0 ){
  A t3(t2.size(),0.0);
  F ckk = 0.0, f1,f2,cc,be;
  int n1 = t1.size();
  for (int k = 1; k <= nn; ++k){
    for (int j = 1; j <= n1; ++j){
      int i1 = k+j-2;
      int i2 = k-j;
      F f1 = 0.0;
      be = (j-1)*delta;
      if (t1[j-1] > 0.0){
        if (i1+1 <= nl){
          f1 = t2[i1+1-1]*exp(-be);
        }
        f2 = 0.0;
        if (i2 >= 0 and i2+1 <= nl){
          f2 = t2[i2+1-1];
        }
        else if (i2 < 0 and 1-i2 <= nl){
          be = -i2*delta;
          f2 = t2[1-i2-1]*exp(-be);
        }
        cc = t1[j-1]*(f1+f2);
        if (j == 1 or j == n1){ cc/=2; }
        t3[k-1] = t3[k-1]+cc;
      }
    }
    t3[k-1] = t3[k-1]*delta;
    if (t3[k-1] <=1e-30){ t3[k-1] = 0.0; }
    cc = t3[k-1];
    be = (k-1)*delta;
    cc += t3[k-1]*exp(-be);
    if ( k == 1 or k == nn){ cc = cc/2; }
    ckk = ckk + cc;
  }
  return t3;
}


*/











