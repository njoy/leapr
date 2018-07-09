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

    f I1C[7] = { 4.5813e-3, 0.0360768, 0.2659732, 1.2067492, 3.0899424, 
                 3.5156229, 1.0 };
    f I2C[7] = { 3.2411e-4, 3.0153e-3, 0.02658733, 0.15084934, 0.51498869, 
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

    f I1C[9] = { 0.00392377, -0.01647633, 0.02635537, -0.02057706, 
      0.00916281, -0.00157565, 0.00225319, 0.01328592, 0.39894228 };

    f I2C[9] = { -0.00420059, 0.01787654, -0.02895312, 0.02282967, 
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

  auto iVec = ranges::view::iota(0,48) | ranges::view::reverse;

  auto InRange = ranges::view::iota(0,50)
                 | ranges::view::transform([](int x){ return x == 1 ? 1.0 : 0.0; });
  

  auto vec1 = ranges::view::zip(iVec,InRange,InRange|ranges::view::slice(1,ranges::end))
            | ranges::view::transform([](auto i){return std::get<0>(i);}) ;
            //| ranges::to_<std::vector<double>>();

    //std::cout << vec1[0] << "    " << vec1[1] << std::endl << std::endl;
    //

  std::vector<double> In ( 50, 0.0 );
  int i = 49;
  In[49] = 0; In[48] = 1;
  
  for ( size_t i = 48; i > 0; --i ){
    In[i-1] = In[i+1] + 2 * (i+1) * In[i] / x;
    if (In[i-1] >= 1.0e10){ 
      for ( size_t j = i; j < In.size(); ++j ){
        In[j-1]=In[j-1]*1.0e-10;
      } 
    }  
  } 

  auto InCutoff = In | ranges::view::transform( [rat=I1/In[0]](auto x){ 
    return cutoff(x * rat); } );

  auto bCutoff = ranges::view::zip( ranges::view::iota(1,51), InCutoff )
               | ranges::view::transform( [expVal, beta_i](auto t){ 
                   int n = std::get<0>(t);
                   if ( std::get<1>(t) == 0.0 ){ return std::make_pair(0.0,0.0); }
                   return std::make_pair( 
                     cutoff( exp(expVal - n*beta_i*0.5) * std::get<1>(t)),
                     cutoff( exp(expVal + n*beta_i*0.5) * std::get<1>(t)) );} );
 //std::cout << std::setprecision(15) << bCutoff <<std::endl;

  auto bplusCutoff  = bCutoff | ranges::view::keys;
  auto bminusCutoff = bCutoff | ranges::view::values;
  /*
  */

  //std::cout << (bminusCutoff|ranges::view::all) << std::endl;
  for ( size_t i = 1; i < In.size(); ++i ){ 
    In[i] = cutoff( In[i] * I1 / In[0] );
  }
  In[0] = cutoff( I1 );

  //std::cout << std::setprecision(15) << (In|ranges::view::all) <<std::endl;
  for ( size_t n = 0; n < In.size(); ++n ){
    if ( In[n] != 0.0 ){
      bplus[n]  = cutoff( exp( expVal - (n+1) * beta_i * 0.5 ) * In[n] );
      bminus[n] = cutoff( exp( expVal + (n+1) * beta_i * 0.5 ) * In[n] );
    }
  }
  //std::cout << std::endl;
  //std::cout << (bplus|ranges::view::all) << std::endl;
  //std::cout << (bminus|ranges::view::all) << std::endl;

  //std::cout << (bminus|ranges::view::all) << std::endl;
  //std::cout << std::endl;
  return std::make_tuple(I0*exp(expVal),bminus|ranges::view::all,bplus|ranges::view::all);
  //return std::make_tuple(I0*exp(expVal),bminusCutoff,bplusCutoff);
}
