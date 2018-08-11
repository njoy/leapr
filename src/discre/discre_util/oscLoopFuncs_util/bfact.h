#include <iostream>
#include <vector>
#include <range/v3/all.hpp>

template <typename F>
F cutoff( F a ){ return a < 1.0e-30 ? 0.0 : a; }

template <typename F, typename A>
auto bfact( const F& x, const F& dwc, const F& beta_i, A& bplus, A& bminus ){
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
    
  F I0, I1, expVal, y = x / 3.75;

  // I_0(x) and I_1(x) are computed first, and then reverse recursion is used
  // to get the rest of the terms (check out numerical recipes for more info).
  // We treat big x and small x differently. Large x is taking advantage of
  // the e^x term, which can be seen when we add x to the var expVal.
  if ( y <= 1.0 ){

    // Compute I0 and I1 via series expansion (manual pg. 714 top)
    F I1C[7] = { 4.5813e-3, 0.0360768, 0.2659732, 1.2067492, 3.0899424, 
                 3.5156229, 1.0 };
    F I2C[7] = { 3.2411e-4, 3.0153e-3, 0.02658733, 0.15084934, 0.51498869, 
                 0.87890594, 0.5};

    auto yVec = ranges::view::iota(0,7) 
              | ranges::view::reverse
              | ranges::view::transform([y](auto num){
                  return std::pow(y,2*num); } );
    auto both = ranges::view::zip(I1C, I2C, yVec) 
              | ranges::view::transform( [](auto t){ 
                  return std::make_pair(std::get<0>(t)*std::get<2>(t),
                         std::get<1>(t)*std::get<2>(t)); } );
    I0 = ranges::accumulate( both|ranges::view::keys,   0.0 );
    I1 = ranges::accumulate( both|ranges::view::values, 0.0 ) * x;

    expVal = -dwc;
  } 
  if ( y > 1.0 ) {

    // Compute I0 and I1 via series expansion (manual pg. 714 top)
    F I1C[9] = { 0.00392377, -0.01647633, 0.02635537, -0.02057706, 
      0.00916281, -0.00157565, 0.00225319, 0.01328592, 0.39894228 };

    F I2C[9] = { -0.00420059, 0.01787654, -0.02895312, 0.02282967, 
      -0.01031555, 0.00163801, -0.00362018, -0.03988024, 0.39894228 };

    auto yVec = ranges::view::iota(0,9) 
              | ranges::view::reverse
              | ranges::view::transform([inv_y=1.0/y](auto num){
                  return std::pow(inv_y,num); } );
    auto both = ranges::view::zip(I1C, I2C, yVec) 
              | ranges::view::transform( [](auto t){ 
                  return std::make_pair(std::get<0>(t)*std::get<2>(t),
                         std::get<1>(t)*std::get<2>(t)); } );
    I0 = ranges::accumulate( both|ranges::view::keys,   0.0 )/sqrt(x);
    I1 = ranges::accumulate( both|ranges::view::values, 0.0 )/sqrt(x);

    expVal = -dwc + x;
  }

  // generate higher orders by reverse recursion
  
  std::vector<double> In ( 50, 0.0 );
  // The "In" vector will be populated with Modified Bessel Functions of the 
  // first kind, all evaluated at x, such that In[0] = I1(x), In[1] = I2(x),
  // etc, where I's subscript denotes the order of the function.
  
  int nmax=50, i = 49;
  In[nmax-1] = 0; In[nmax-2] = 1;
  while (i > 1){
    In[i-2] = In[i] + 2 * i * In[i-1] / x;
    i = i - 1;
    if (In[i-1] >= 1.0e10){ 
      for ( auto j = i; j < nmax; ++j ){
        In[j-1]=In[j-1]/1.0e10;
      } 
    }  
  } 

  for ( auto i = 1; i < nmax; ++i ){ 
    In[i] = cutoff( In[i] * I1 / In[0] );
  }
  In[0] = cutoff( I1 );
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

  /* blus and bminus are being populated with the A_in values for Eq. 537.
   * The sum is for n ranging from -infty to infty, but we approximate by 
   * having n range between 0 and 50, and putting +/- vals in bplus/minus.
   */

  for ( auto n = 0; n < nmax; ++n ){
    bplus[n]  = cutoff( exp( expVal - (n+1) * beta_i / 2 ) * In[n] );
    bminus[n] = cutoff( exp( expVal + (n+1) * beta_i / 2 ) * In[n] );
  }

  // Return bzero value
  return I0 * exp( expVal );
}

