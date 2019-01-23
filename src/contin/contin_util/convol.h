#include <iostream>
#include <range/v3/all.hpp>


template <typename F, typename A>
auto getTnVal( int j, const A& tn, const F& delta ){
  if (abs(j) >= int(tn.size())){ return 0.0; }
  if (j < 0){ return tn[abs(j)]*exp(j*delta); }
  return tn[abs(j)];
}

template <typename F, typename A>
auto getConvolAtPoint( int i, F delta, const A& t1, const A& t2){
  int len_t1 = int(t1.size());
  F sumVal = 0.0;
  for ( int j = -len_t1+1; j < len_t1; ++j ){
    if ( i - j  >= int(t2.size()) ){ return sumVal; }
    F value = getTnVal( j, t1, delta ) * getTnVal( i-j, t2, delta );
    sumVal += (abs(j) == len_t1-1) ? value * 0.5 * delta : value * delta;
  }
  return sumVal;
}

template <typename F, typename A>
auto convol( const A& t1, const A& t2, const F& delta, int len_t3){
  A t3(t2.size(),0.0);
  for ( int i = 0; i < len_t3; ++i ){    
      t3[i] = getConvolAtPoint(i,delta,t1,t2);
  } 
  return t3;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


template <typename F, typename A>
auto getTnVal2( int j, const A& tn, const F& delta ){
  if (abs(j) >= int(tn.size())){ return 0.0; }
  if (j < 0){ return tn[abs(j)]*exp(j*delta); }
  return tn[abs(j)];
}

template <typename F, typename A>
auto getConvolAtPoint2( int i, F delta, const A& t1, const A& t2){
  int len_t1 = int(t1.size());
  F sumVal = 0.0;
  for ( int j = -len_t1+1; j < len_t1; ++j ){
    if ( i - j  >= int(t2.size()) ){ return sumVal; }
    F value = getTnVal2( j, t1, delta ) * getTnVal2( i-j, t2, delta );
    sumVal += (abs(j) == len_t1-1) ? value * 0.5 * delta : value * delta;
  }
  return sumVal;
}


template <typename F, typename A>
auto convol2( const A& t1, const A& t2, const F& delta, int len_t3){
  A t3(t2.size(),0.0);
  for ( int i = 0; i < len_t3; ++i ){    
      t3[i] = getConvolAtPoint2(i,delta,t1,t2);
  } 
  return t3;
}







