//#include <iostream>
//#include <stdlib.h> // abs()
#include <iostream>
#include <range/v3/all.hpp>

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

template <typename F, typename A>
auto getConvol2AtPoint( int i, int nl, F delta, const A& t1, const A& t2, int len_t1 ){
    F sumVal = 0.0;
    for ( int j = -len_t1+1; j < len_t1; ++j ){  // j iterates through t1, t2
      F val = 1.0;
      if ( i - j >= int(t2.size()) ){ return sumVal; }
      if (j   < 0){ val *= exp(j*delta);}
      if (i-j < 0){ val *= exp((i-j)*delta); }

      F toAdd = t1[abs(j)]*t2[abs(i-j)]*val;
      //std::cout << "-------------  " <<i-j << "    " << toAdd << std::endl;
      sumVal += ( j == -len_t1+1 or j == len_t1-1 ) ? 0.5*toAdd : toAdd;

    } // for
    return sumVal;
    if ( nl+i> 0){sumVal += + t1[0] + t2[1]+delta;}
}



template <typename F, typename A>
auto convol( const A& t1, const A& t2, const F& delta, const int& nl=0,
  const int& n1=0, const int& nn=0 ){

  A t3(t2.size(),0.0);

  int len_t1 = int(t1.size()), len_t2 = int(t2.size());
  F f1, f2;
  len_t2 = nn;
  len_t1 = n1;

//  std::cout << getConvolAtPoint(0,nl,delta,t1,t2,len_t1)-getConvol2AtPoint(0,nl,delta,t1,t2,len_t1) << std::endl;
  /*
  */

  for ( int i = 0; i < len_t2; ++i ){    // i iterates through t3
      t3[i] = getConvolAtPoint(i,nl,delta,t1,t2,len_t1);
      //std::cout << i << std::endl;
      //std::cout << i << "     " << t3[i] << "    ";
      t3[i] = getConvol2AtPoint(i,nl,delta,t1,t2,len_t1);
      //std::cout << t3[i] << std::endl;
      t3[i] = ( t3[i] * delta < 1e-30 ) ? 0 : t3[i] * delta;

  } // for
  return t3;

}


/*
template <typename F, typename A>
auto getTnVal( int j, const A& tn, const F& delta ){
  if (abs(j) >= int(tn.size())){ return 0.0; }
  if (j < 0){ return tn[abs(j)]*exp(j*delta); }
  return tn[abs(j)];
}


template <typename A, typename Range, typename B>
auto getConvolAtPoint( int i, const A& t1, const B& t2, const Range& betaGrid ){
  int len_t1 = int(t1.size());
  double sumVal = 0.0, val_L, val_R;
  for (int j = -len_t1+1; j < len_t1-1; ++j){
    if (i-j  >= int(t2.size())){ return sumVal; }
    auto delta = abs(betaGrid[abs(j+1)] - betaGrid[abs(j)]);
    val_L = getTnVal(j,  t1,delta)*getTnVal(i-j,    t2,delta);
    val_R = getTnVal(j+1,t1,delta)*getTnVal(i-(j+1),t2,delta);
    sumVal += (val_L + val_R)*0.5*delta;
  }
  return sumVal;
}


template <typename Range, typename Array>
auto convol( const Range& t1, const Array& t2, int len_t3, const Array& betaGrid){
  Array t3(t2.size(),0.0);
 
  for (int i = 0; i<len_t3; ++i){
    t3[i] = getConvolAtPoint(i,t1,t2,betaGrid);
  }
  return t3;
}

*/



/*
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
*/


