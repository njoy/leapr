#include <iostream>
#include <range/v3/all.hpp>

template <typename F, typename A>
auto getConvolAtPoint( int i, F delta, const A& t1, const A& t2, int len_t1 ){
  /* t3 = t1 * t2
   * t3(i) = sum_j t1(j) t2(i-j)
   *
   * also t(x) = exp(-x)*t(-x)
   *
   * t3(1) = t1(-4) t2( 5)  -->  t1(4)*exp(-4delta) t2(5) --> t1(4)*t2(4+1)*exp(-4delta)
   *       + t1(-3) t2( 4)  -->  t1(3)*exp(-3delta) t2(4) --> t1(3)*t2(3+1)*exp(-3delta)
   *       + t1(-2) t2( 3)  -->  t1(2)*exp(-2delta) t2(3) --> t1(2)*t2(2+1)*exp(-2delta)
   *       + t1(-3) t2( 2)  -->  t1(1)*exp(-1delta) t2(2) --> t1(1)*t2(1+1)*exp(-1delta)
   *       + t1( 0) t2( 1)  -->  t1(0)*exp(-0delta) t2(1) --> t1(0)*t2(0+1)*exp(-0delta)
   *       + t1( 1) t2( 0)  -->  t1(1) t2(0)*exp(-0delta) --> t1(1)*t2(|-1+1|)*exp(-1delta)
   *       + t1( 2) t2(-1)  -->  t1(2) t2(1)*exp(-1delta) --> t1(2)*t2(|-2+1|)*exp(-2delta)
   *       + t1( 3) t2(-2)  -->  t1(3) t2(2)*exp(-2delta) --> t1(3)*t2(|-3+1|)*exp(-3delta)
   *       + t1( 4) t2(-3)  -->  t1(4) t2(3)*exp(-3delta) --> t1(4)*t2(|-4+1|)*exp(-4delta)
   */
  F sumVal = 0.0;
  for ( int j = -len_t1+1; j < len_t1; ++j ){  // j iterates through t1, t2
    if ( i - j >= int(t2.size()) ){ return sumVal; }
    F val = 1.0;
    if (j   < 0){ val *= exp(j*delta);}
    if (i-j < 0){ val *= exp((i-j)*delta); }

    F toAdd = t1[abs(j)]*t2[abs(i-j)]*val;
    sumVal += ( j == -len_t1+1 or j == len_t1-1 ) ? 0.5*toAdd : toAdd;

  } // for
  return sumVal;
}



template <typename F, typename A>
auto convol( const A& t1, const A& t2, const F& delta, const int& nl=0,
  const int& n1=0, const int& nn=0 ){

  A t3(t2.size(),0.0);

  int len_t1 = int(t1.size()), len_t2 = int(t2.size());
  F f1, f2;
  len_t2 = nn;
  len_t1 = n1;

  for ( int i = 0; i < len_t2; ++i ){    // i iterates through t3
      t3[i] = getConvolAtPoint(i,delta,t1,t2,len_t1);
      t3[i] = ( t3[i] * delta < 1e-30 ) ? 0 : t3[i] * delta;

  } // for
  return t3;

  std::cout << nl << std::endl;
}

