#include <iostream> 
#include <vector>
#include <cmath>

auto bessel_K1_gen(double x ){
  double v,u,bi1,bi3;
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

  double besk1;
  
  if (x <= 1 ){
    v=c0*x;
    u=v*v;
    bi1=(((((((((c1*u+c2)*u+c3)*u+c4)*u+c5)*u+c6)*u
                     +c7)*u+c8)*u+c9)*u+c10)*v;
    bi3=(((((((((c11*u+c12)*u+c13)*u+c14)*u+c15)*u
                     +c16)*u+c17)*u+c18)*u+c19)*u+c20);
    besk1=c21/x+ bi1*(log(c22*x)+c23)-v*bi3;
  }
  else {
    u=c24/x;
    bi3=((((((((((((-c25*u+c26)*u-c27)*u+c28)*u-c29)*u+c30)*u-c31)*u
                          +c32)*u-c33)*u+c34)*u-c35)*u+c36)*u+c37);
    besk1=sqrt(u)*bi3;
  }
  return besk1;
}





