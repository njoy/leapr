#include <iostream> 
#include <vector>
#include <cmath>

auto bessel_K1_gen(double x ){
  double c0  = 0.125;
  double c1  = 0.442850424;
  double c2  = 0.584115288;
  double c3  = 6.070134559;
  double c4  = 17.864913364;
  double c5  = 48.858995315;
  double c6  = 90.924600045;
  double c7  = 113.795967431;
  double c8  = 85.331474517;
  double c9  = 32.00008698;
  double c10 = 3.999998802;
  double c11 = 1.304923514;
  double c12 = 1.47785657;
  double c13 = 16.402802501;
  double c14 = 44.732901977;
  double c15 = 115.837493464;
  double c16 = 198.437197312;
  double c17 = 222.869709703;
  double c18 = 142.216613971;
  double c19 = 40.000262262;
  double c20 = 1.999996391;
  double c21 = 1.;
  double c22 = 0.5;
  double c23 = 0.5772156649;
  double c24 = 1.;
  double c25 = 0.0108241775;
  double c26 = 0.0788000118;
  double c27 = 0.2581303765;
  double c28 = 0.5050238576;
  double c29 = 0.663229543;
  double c30 = 0.6283380681;
  double c31 = 0.4594342117;
  double c32 = 0.2847618149;
  double c33 = 0.1736431637;
  double c34 = 0.1280426636;
  double c35 = 0.1468582957;
  double c36 = 0.4699927013;
  double c37 = 1.2533141373;
  double v,u,bi1,bi3;

  std::vector<double> constVec; 
  if (x <= 1 ){
    v = c0 * x;

    bi1 = c1*v*v;
    constVec = {c2,c3,c4,c5,c6,c7,c8,c9};
    for ( auto& entry : constVec ){ bi1 = (bi1 + entry)*v*v; }
    bi1 = ( bi1 + c10 ) * v;

    bi3 = c11*v*v;
    constVec = {c12,c13,c14,c15,c16,c17,c18,c19};
    for ( auto& entry : constVec ){ bi3 = (bi3 + entry)*v*v; }
    bi3 = bi3 + c20;
    
    return c21 / x + bi1 * ( log( c22 * x ) + c23 ) - v * bi3;

  }
  else {
    u = c24 / x;
    bi3 = -c25*u;
    int i = 1;
    constVec = {c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36};
    for (auto& entry : constVec ){
      bi3 = ( bi3 + i*entry ) * u;
      i = -i;
    }
    return sqrt( u ) * ( bi3 + c37 );
  }
}





