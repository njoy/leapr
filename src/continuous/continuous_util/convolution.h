#include <cmath>

template <typename Float, typename Range>
auto getConvolAtPoint( int i, const Float& delta, const Range& t1, 
  const Range& t2 ){
  Float sumVal = 0.0, expVal;
  int len_t1 = t1.size();
  for ( int j = -len_t1+1; j < len_t1; ++j ){  // j iterates through t1, t2
    if ( i - j >= int(t2.size()) ){ return sumVal; }
    expVal = 1.0;
    if (j < 0){ expVal = exp(j*delta); }
    if (i < j){ expVal = exp((i-j)*delta); }

    sumVal += ( j == -len_t1+1 or j == len_t1-1 ) ? 
                0.5*(t1[std::fabs(j)]*t2[std::fabs(i-j)]*expVal)
              :     (t1[std::fabs(j)]*t2[std::fabs(i-j)]*expVal);
  } 
  return sumVal;
}


template <typename Float, typename Range>
auto convolution( const Range& t1, const Range& t2, const Float& delta, const int nn){
  Range t3(nn,0.0);
  for ( int i = 0; i < nn; ++i ){   
    t3[i] = getConvolAtPoint(i,delta,t1,t2) * delta;
  } 
  return t3;
}

