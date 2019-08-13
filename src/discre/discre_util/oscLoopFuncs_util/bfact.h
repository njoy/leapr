#include <range/v3/all.hpp>

template <typename Float>
auto cutoff( Float a ){ return a < 1.0e-30 ? 0.0 : a; }

template <typename Float, typename Range>
auto bfact( const Float& x, const Float& dwc, const Float& beta_i, 
  Range& bplus, Range& bminus ){
  using std::exp;
   
  Float I0, I1, expVal = -dwc, y = x / 3.75;

  if ( y <= 1.0 ){

    Range 
      I1Consts = {0.0360768, 0.2659732, 1.2067492, 3.0899424, 3.5156229, 1.0},
      I2Consts = {0.0030153, 0.02658733,0.15084934,0.51498869,0.87890594,0.5};
    I0 = 0.0045813;
    I1 = 0.00032411;
    for ( auto entry : I1Consts ){ I0 = I0 * y * y + entry; }
    for ( auto entry : I2Consts ){ I1 = I1 * y * y + entry; }

    I1 = I1 * x; 

  } 
  if ( y > 1.0 ) {

    // Compute I0 and I1 via series expansion (manual pg. 714 top)
    std::vector<double> I1Consts { -0.01647633, 0.02635537, -0.02057706, 
      0.00916281, -0.00157565, 0.00225319, 0.01328592, 0.39894228 };
    std::vector<double> I2Consts { 0.01787654, -0.02895312, 0.02282967, 
      -0.01031555, 0.00163801, -0.00362018, -0.03988024, 0.39894228 };
    I0 = 0.00392377;
    I1 = -0.00420059;

    for ( auto entry : I1Consts){ I0 = I0 / y + entry; }
    for ( auto entry : I2Consts){ I1 = I1 / y + entry; }

    I0 /= sqrt(x);
    I1 /= sqrt(x);

    expVal += x;
  }

  Range In ( 50, 0.0 );
  
  int nmax=50, i = 49;
  In[nmax-1] = 0; In[nmax-2] = 1;
  while (i > 1){
    In[i-2] = In[i] + 2 * i * In[i-1] / x;
    i = i - 1;
    if (In[i-1] >= 1.0e10){ 
      for (auto j = i-1; j < nmax-1; ++j){ In[j] *= 1e-10; } 
    }  
  } 

  Float fraction = I1/In[0];
  for ( auto i = 0; i < nmax; ++i ){ 
    In[i] = cutoff( In[i] * fraction ); 
  }

  for ( auto n = 0; n < nmax; ++n ){
    bplus[n]  = cutoff( exp(expVal - (n+1) * beta_i * 0.5) * In[n] );
    bminus[n] = cutoff( exp(expVal + (n+1) * beta_i * 0.5) * In[n] );
  }

  return I0 * exp(expVal);
}
