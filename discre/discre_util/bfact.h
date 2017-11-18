#include <iostream>
#include <vector>

auto bfact( const double& x, const double& dwc, const double& beta_i, 
  std::vector<double>& bplus, std::vector<double>& bminus ){
    
   // compute bessi0
  double bessi0;
  double y = x / 3.75;
  if ( y < 1.0 ){
    double u = y * y;
    std::vector<double> constVec { 0.0360768, 0.2659732, 1.2067492, 3.0899424, 
      3.5156229, 1.0 };
    bessi0 = 0.0045813;
    for ( auto entry : constVec ){ bessi0 = bessi0 * u + entry; }
  } else {
    double v = 1 / y;
    std::vector<double> constVec { -0.01647633, 0.02635537, -0.02057706, 
      0.00916281, -0.00157565, 0.00225319, 0.01328592, 0.39894228 };
    bessi0 = 0.00392377;
    for ( auto entry : constVec ){ bessi0 = bessi0 * v + entry; }
    bessi0 = bessi0/sqrt(x);
  }

  // compute bessi1
  double bessi1;
  if ( y <= 1.0 ){
    double u = y * y;
    std::vector<double> constVec { 0.00301532, 0.02658733, 0.15084934, 
      0.51498869, 0.87890594, 0.5 };
    bessi1 = 0.00032411;
    for ( auto entry : constVec ){ bessi1 = bessi1 * u + entry; }
    bessi1 = bessi1 * x; 
  } else {
    double v = 1 / y;
    std::vector<double> constVec { 0.01787654, -0.02895312, 0.02282967, 
      -0.01031555, 0.00163801, -0.00362018, -0.03988024, 0.39894228 };
    bessi1 = -0.00420059;
    for ( auto entry : constVec ){ bessi1 = bessi1 * v + entry; }
    bessi1 = bessi1 / sqrt(x);
  }

  // generate higher orders by reverse recursion
  std::vector<double> bn ( 50, 0.0 );
  int imax=50;
  bn[imax-1] = 0;
  bn[imax-2] = 1;
  int i = imax - 1;
  while (i > 0){
    bn[i-2]=bn[i]+i*(2/x)*bn[i-1];
    i = i - 1;
    if (bn[i-1] >= 1e10){ 
      for ( auto j = i; j < imax; ++j ){
        bn[j-1]=bn[j-1]/1.0e10;
      } 
    }  
  } 

  auto rat = bessi1/bn[0];
  for ( auto i = 0; i < imax; ++i ){ 
    bn[i] = bn[i] * rat < 1.0e-30 ? 0.0 : bn[i] * rat;
  }

  // apply exponential terms to bessel functions
  double bzero, arg, expVal;
  expVal = y <= 1.0 ? -dwc : -dwc + x;
  bzero=bessi0*exp(expVal);
  for ( auto i = 0; i < imax; ++i ){
    bplus[i]  = exp( expVal - (i+1)*beta_i/2 ) * bn[i];
    bminus[i] = exp( expVal + (i+1)*beta_i/2 ) * bn[i];
    if (bminus[i] < 1e-30 ) { bminus[i] = 0; }
    if (bplus[i]  < 1e-30 ) { bplus[i]  = 0; }
  }
  return bzero;
}
