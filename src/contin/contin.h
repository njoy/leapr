#include "contin/contin_util/start.h"
#include "contin/contin_util/convol.h"
#include "generalTools/interpolate.h"
#include <range/v3/all.hpp>

template <typename Range, typename Float>
auto contin(int nphon, const Float& delta, const Float& continWgt, 
            const Range& rho, const Range& alpha, const Range& beta, Range& sab ){
  using std::exp;

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


  auto lambda_alpha     = ranges::view::iota(0,int(alpha.size()))
                        | ranges::view::transform([&](auto i){ 
                            return lambda_s*alpha[i]; });
  auto exp_lambda_alpha = lambda_alpha 
                        | ranges::view::transform([](auto lambda_alpha){
                            return exp(-lambda_alpha); });

  for( int n = 0; n < nphon; ++n ){
    if ( n > 0 ){ 
      npn = t1.size()+npl-1;
      tnow = convol(t1, tlast, delta, npn); 
    }
    double inv_n = 1.0/(n+1);

    for( size_t a = 0; a < alpha.size(); ++a ){
      xa[a] *= lambda_alpha[a] * inv_n;
      exx    = exp_lambda_alpha[a]*xa[a];
 
      for( size_t b = 0; b < beta.size(); ++b ){
        add = exx * interpolate(tnow, beta[b], delta);
        sab[b+a*beta.size()] += add < 1e-30 ? 0 : add;
      } 
    } 

    if ( n == 0 ){ continue; }
    for( size_t i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }
    npl = npn;
  } 

  return std::make_tuple(lambda_s,t_bar);
}



