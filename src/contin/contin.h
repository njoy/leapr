#include "contin/contin_util/start.h"
#include "contin/contin_util/convol.h"
#include "generalTools/interpolate.h"
#include <range/v3/all.hpp>

template <typename Range, typename Float>
auto contin(int nphon, Float& delta, const Float& continWgt, 
  const Float& scaling, const Float& tev, const Float& sc, Range rho, 
  Range& alpha, Range& beta, Range& symSab ){
  using std::exp;

  for ( auto& a : alpha ){ a *= scaling; }
  for ( auto& b : beta  ){ b *= sc;      }

  delta /= tev;

  Range betaGrid = ranges::view::iota(0,int(rho.size())) 
                 | ranges::view::transform([delta](auto x){return x*delta;});
    
  auto lambda_s_t_bar = start( rho, continWgt, betaGrid );
  Float lambda_s = std::get<0>(lambda_s_t_bar),
        t_bar    = std::get<1>(lambda_s_t_bar);
  Range t1       = std::get<2>(lambda_s_t_bar);
  
  Range xa(alpha.size(),1.0), tnow(nphon*t1.size(),0.0), tlast(nphon*t1.size(),0.0);
  std::copy( t1.begin(), t1.begin() + t1.size(), tlast.begin() );
  std::copy( t1.begin(), t1.begin() + t1.size(), tnow.begin()  );

  Float add, exx;
  
  size_t npn = t1.size(), npl = t1.size();

  Range exp_lambda_alpha = ranges::view::iota(0,int(alpha.size()))
                         | ranges::view::transform([&](auto i){ 
                             return exp(-lambda_s*alpha[i]);
                           });

  for( int n = 0; n < nphon; ++n ){
    if ( n > 0 ){ 
      npn = t1.size()+npl-1;
      tnow = convol2(t1, tlast, delta, npl, npn); 
    }

    for( size_t a = 0; a < alpha.size(); ++a ){
      xa[a] *= lambda_s * alpha[a] / ( n + 1 );
      exx    = exp_lambda_alpha[a]*xa[a];
      Range betaGrid2 = ranges::view::iota(0,int(tnow.size())) 
                      | ranges::view::transform([delta](auto x){return x*delta;});
 
      for( size_t b = 0; b < beta.size(); ++b ){
        add = exx * interpolate(tnow, beta[b], betaGrid2);
        symSab[b+a*beta.size()] += add < 1e-30 ? 0 : add;
      } 
    } 

    if ( n == 0 ){ continue; }
    for( size_t i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }
    npl = npn;
  } 

  return std::make_tuple(lambda_s,t_bar);
}



