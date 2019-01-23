#include <iostream>
#include <range/v3/all.hpp>


template <typename F, typename A>
auto getTnVal( int i, const A& tn, const F& delta ){
  if (abs(i) >= int(tn.size())){ return 0.0; }
  if (i < 0){ return tn[abs(i)]*exp(i*delta); }
  return tn[abs(i)];

}
 

template <typename F, typename A>
auto getConvolAtPointg int i, F delta, const A& t1, const A& t2){
  int len_t1 = int(t1.size());
  F sumVal = 0.0;
  for ( int j = -len_t1+1; j < len_t1; ++j ){
    if ( i - j  >= int(t2.size()) ){ return sumVal; }
    F value = getTnVal( j, t1, delta ) * getTnVal( i-j, t2, delta );
    sumVal += (abs(j) == len_t1-1) ? value * 0.5 : value;
  }
  return sumVal;
}




template <typename F, typename A>
auto convol( const A& t1, const A& t2, const F& delta, int len_t3){
  A t3(t2.size(),0.0);
  F f1, f2;
  for ( int i = 0; i < len_t3; ++i ){    
      t3[i] = getConvolAtPoint(i,delta,t1,t2)*delta;
  } 
  return t3;
}






