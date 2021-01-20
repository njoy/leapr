#include "continuous/continuous_util/start.h"
#include "continuous/continuous_util/convolution.h"
#include "continuous/continuous_util/checkMoments.h"
#include "generalTools/tools.h"
#include <range/v3/all.hpp>


template <typename Range, typename Float>
auto continuous(int nphon, const Float& delta, const Float& continWgt, 
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

  Float add, exx, inv_n;
  
  size_t nNext = t1.size(), nLast = t1.size();

  auto lambda_alpha     = ranges::view::iota(0,int(alpha.size()))
                        | ranges::view::transform([&](auto i){ 
                            return lambda_s*alpha[i]; });
  auto exp_lambda_alpha = lambda_alpha 
                        | ranges::view::transform([](auto lambda_alpha){
                            return exp(-lambda_alpha); });

  std::vector<int> maxt(beta.size(),alpha.size()+1);

  for( int n = 0; n < nphon; ++n ){
    if ( n > 0 ){ 
      nNext = t1.size()+nLast-1;
      tnow = convolution(t1, tlast, delta, nNext); 
    }
    inv_n = 1.0/(n+1);

    for( size_t a = 0; a < alpha.size(); ++a ){
      xa[a] *= lambda_alpha[a] * inv_n;
      exx    = exp_lambda_alpha[a]*xa[a];
 
      for( size_t b = 0; b < beta.size(); ++b ){
        add = exx * interpolate(tnow, beta[b], delta, nNext);
        if (add < 1e-30){ add = 0; }
        sab[b+a*beta.size()] += add;
        if (n == nphon - 1 and sab[b+a*beta.size()] > 0 and 
             int(a+1) < maxt[b] and 1000*add > sab[b+a*beta.size()] ){
          maxt[b] = a+1; 
        }
      } 
    } 

    if ( n == 0 ){ continue; }
    for( size_t i = 0; i < nNext; ++i ){ tlast[i] = tnow[i]; }
    nLast = nNext;
  } 

  checkMoments(alpha,beta,maxt,lambda_s,continWgt,t_bar,sab);

  return std::make_tuple(lambda_s,t_bar);
}



