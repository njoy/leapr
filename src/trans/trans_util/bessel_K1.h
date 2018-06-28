#ifndef TRANS_TRANSUTIL_BESSELK1GEN
#define TRANS_TRANSUTIL_BESSELK1GEN
#include <vector>
#include <iostream> 
#include <cmath>

template <typename floatT>
inline floatT bessel_K1_gen(floatT x ){
  /* Computes the modified bessel function of the second kind, K1(x). For 
   * x values greater than 1, the exponential term is omitted.
   */ 
  floatT c0 = 0.125,         c1 = 0.442850424,   c2  = 0.584115288,
         c3 = 6.070134559,   c4 = 17.864913364,  c5  = 48.858995315,
         c6 = 90.924600045,  c7 = 113.795967431, c8  = 85.331474517,
         c9 = 32.00008698,  c10 = 3.999998802,   c11 = 1.304923514,
        c12 = 1.47785657,   c13 = 16.402802501,  c14 = 44.732901977,
        c15 = 115.83749346, c16 = 198.437197312, c17 = 222.869709703,
        c18 = 142.216613971,c19 = 40.000262262,  c20 = 1.999996391,
        c21 = 1.,           c22 = 0.5,           c23 = 0.5772156649,
        c24 = 1.,           c25 = 0.0108241775,  c26 = 0.0788000118,
        c27 = 0.2581303765, c28 = 0.5050238576,  c29 = 0.663229543,
        c30 = 0.6283380681, c31 = 0.4594342117,  c32 = 0.2847618149,
        c33 = 0.1736431637, c34 = 0.1280426636,  c35 = 0.1468582957,
        c36 = 0.4699927013, c37 = 1.2533141373;
  floatT v,bi1,bi3;

  std::vector<floatT> constVec; 
  if (x <= 1 ){
    v = c0 * x;

    bi1 = c1 * v * v;
    constVec = {c2,c3,c4,c5,c6,c7,c8,c9};
    for ( const auto& entry : constVec ){ bi1 = (bi1 + entry)*v*v; }
    bi1 = ( bi1 + c10 ) * v;

    bi3 = c11 * v * v;
    constVec = {c12,c13,c14,c15,c16,c17,c18,c19};
    for ( const auto& entry : constVec ){ bi3 = (bi3 + entry)*v*v; }
    bi3 = bi3 + c20;
    
    return c21 / x + bi1 * ( log( c22 * x ) + c23 ) - v * bi3;

  } else {
    v = c24 / x;
    bi3 = -c25 * v;
    int i = 1;
    constVec = {c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36};
    for (const auto& entry : constVec ){
      bi3 = ( bi3 + i*entry ) * v;
      i = -i;
    }
    return sqrt( v ) * ( bi3 + c37 );
  }
}

#endif



