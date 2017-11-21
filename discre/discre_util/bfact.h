#include <iostream>
#include <vector>

auto cutoff( double a ){
  return a < 1.0e-30 ? 0.0 : a;
}

auto bfact( const double& x, const double& dwc, const double& beta_i, 
  std::vector<double>& bplus, std::vector<double>& bminus ){
  /* bfact is meant to assist with the evaluation of Eq. 537, by evaluating
   * the bessel function In and all exponential terms. 
   *
   * This calculates the bessel function terms for the discrete oscillator 
   * terms. The bessel functions I_n are evaluated at point x and put into
   * vector In such that In[0] = I1(x) (modified bessel function of first kind
   * of order 1), In[1] = I2(x) (modified bessel function of first kind of
   * order 2), etc. 
   *
   * The x that is input here is equal to 
   *             alpha * weight / ( beta * sinh( beta / 2 ) )
   * to match Eq. 537. 
   *
   * The dwc vector is populated with entries of dwc[i] = alpha * lambda_i,
   * where lambda_i was calculated in prepareParams.h according to Eq. 538.
   */
    
  double bessi0;
  double bessi1, expVal;
  double y = x / 3.75;
   // compute bessi0
  // compute bessi1
  if ( y <= 1.0 ){

    // Compute bessi1 and bessi1
    std::vector<double> I1Consts = { 0.0360768, 0.2659732, 1.2067492, 
      3.0899424, 3.5156229, 1.0 };
    std::vector<double> I2Consts = {0.0030153, 0.02658733, 0.15084934, 
      0.51498869, 0.87890594, 0.5};
    bessi0 = 0.0045813;
    bessi1 = 0.00032411;
    for ( auto entry : I1Consts ){ bessi0 = bessi0 * y * y + entry; }
    for ( auto entry : I2Consts ){ bessi1 = bessi1 * y * y + entry; }

    bessi1 = bessi1 * x; 

    expVal = -dwc;
  } 
  if ( y > 1.0 ) {

    // Compute bessi0 and bessi1
    std::vector<double> I1Consts { -0.01647633, 0.02635537, -0.02057706, 
      0.00916281, -0.00157565, 0.00225319, 0.01328592, 0.39894228 };
    std::vector<double> I2Consts { 0.01787654, -0.02895312, 0.02282967, 
      -0.01031555, 0.00163801, -0.00362018, -0.03988024, 0.39894228 };
    bessi0 = 0.00392377;
    bessi1 = -0.00420059;

    for ( auto entry : I1Consts){ bessi0 = bessi0 / y + entry; }
    for ( auto entry : I2Consts){ bessi1 = bessi1 / y + entry; }

    bessi0 /= sqrt(x);
    bessi1 /= sqrt(x);

    expVal = -dwc + x;
  }

  // generate higher orders by reverse recursion
  
  std::vector<double> In ( 50, 0.0 );
  // The "In" vector will be populated with Modified Bessel Functions of the 
  // first kind, all evaluated at x, such that In[0] = I1(x), In[1] = I2(x),
  // etc, where I's subscript denotes the order of the function.
  
  int imax=50, i = 49;
  In[imax-1] = 0; In[imax-2] = 1;
  while (i > 0){
    In[i-2] = In[i] + 2 * i * In[i-1] / x;
    i = i - 1;
    if (In[i-1] >= 1.0e10){ 
      for ( auto j = i; j < imax; ++j ){
        In[j-1]=In[j-1]/1.0e10;
      } 
    }  
  } 

  for ( auto i = 1; i < imax; ++i ){ 
    In[i] = cutoff( In[i] * bessi1 / In[0] );
  }
  In[0] = cutoff( bessi1 );
  // bessel function vector In is finished being populated


  /* Having calculated In(alpha*weight/(beta*sinh(beta/2))) = In(x), we apply
   * exponential terms to it. expVal = -dwc or -dwc + x. 
   * dwc = alpha * lambda_i, where lambda_i was calculated in prepareParams.h
   * according to Eq. 538. 
   * Additionally, the exp( -n * beta / 2 ) term in Eq. 537 is evaluated.
   * In total, bplus and bminus will be populated with 
   *      exp(-a*lambda_i) * In( a*weight/(b*sinh(b/2)) ) * exp(-n*b/2)
   * for positive and negative values of n, respectively. (b = beta, a = alpha)
   */
  for ( auto i = 0; i < imax; ++i ){
    bplus[i]  = cutoff( exp( expVal - (i+1) * beta_i / 2 ) * In[i] );
    bminus[i] = cutoff( exp( expVal + (i+1) * beta_i / 2 ) * In[i] );
  }
  // Return bzero value
  return bessi0 * exp( expVal );
}
