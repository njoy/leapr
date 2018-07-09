#include <iostream>
#include <vector>
#include <range/v3/all.hpp>

template <typename floatT>
floatT cutoff( floatT a ){ return a < 1.0e-30 ? 0.0 : a; }

template <typename f, typename arrayT>
auto bfact( const f& x, const f& dwc, const f& beta_i, 
  arrayT& bplus, arrayT& bminus ){
  f I0;
  f I1, expVal;
  f y = x / 3.75;

  if ( y <= 1.0 ){

    f I1C[6] = { 0.0360768, 0.2659732, 1.2067492, 3.0899424, 3.5156229, 1.0 };
    f I2C[6] = { 3.0153e-3, 0.02658733, 0.15084934, 0.51498869, 0.87890594, 0.5};

    auto I0Vector = ranges::view::concat(ranges::view::single(0.0045813),I1C);
    auto I1Vector = ranges::view::concat(ranges::view::single(0.00032411),I2C);

    auto yVec = ranges::view::iota(0,7) 
              | ranges::view::transform([size=(sizeof(I1C)/sizeof(*I1C))]
                  (auto i){ return 2*(size - i); })
              | ranges::view::transform([y](auto num){
                  return std::pow(y,num); } );
    auto both = ranges::view::zip(I0Vector, I1Vector, yVec) 
              | ranges::view::transform( [](auto t){ 
                  return std::make_pair(std::get<0>(t)*std::get<2>(t),
                         std::get<1>(t)*std::get<2>(t)); } );
    auto I02_1 = ranges::accumulate( both|ranges::view::keys, 0.0 );
    auto I12_1 = ranges::accumulate( both|ranges::view::values, 0.0 );

    std::cout << I02_1 << std::endl;
    std::cout << I12_1 << std::endl;
    std::cout << std::endl;
    



    std::vector<double> I1Consts = { 0.0360768, 0.2659732, 1.2067492, 
      3.0899424, 3.5156229, 1.0 };
    std::vector<double> I2Consts = {0.0030153, 0.02658733, 0.15084934, 
      0.51498869, 0.87890594, 0.5};
    I0 = 0.0045813;
    I1 = 0.00032411;
    for ( auto entry : I1Consts ){ 
      I0 = I0 * y * y + entry; 
    }
    for ( auto entry : I2Consts ){ I1 = I1 * y * y + entry; }
    std::cout << I0 << std::endl;
    std::cout << I1 << std::endl;

    I1 = I1 * x; 

    expVal = -dwc;
  } 
  if ( y > 1.0 ) {

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

    expVal = -dwc + x;
  }

  std::vector<double> In ( 50, 0.0 );
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

  for ( auto n = 0; n < nmax; ++n ){
    bplus[n]  = cutoff( exp( expVal - (n+1) * beta_i / 2 ) * In[n] );
    bminus[n] = cutoff( exp( expVal + (n+1) * beta_i / 2 ) * In[n] );
  }

  return I0 * exp( expVal );
}
