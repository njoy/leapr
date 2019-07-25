#include "contin/contin_util/start.h"
#include "contin/contin_util/convol.h"
#include "generalTools/interpolate.h"
#include <range/v3/all.hpp>

template <typename A, typename F>
auto contin(int nphon, F& delta, const F& tbeta, const F& scaling, const F& tev, 
  const F& sc, A rho, const A& alpha, const A& beta, A& symSab, A betaGrid ){

  for ( F& x : betaGrid ){ x /= tev; }
    
  auto lambda_s_t_eff = start( rho, tbeta, betaGrid );
  F lambda_s = std::get<0>(lambda_s_t_eff),
    t_eff    = std::get<1>(lambda_s_t_eff);
  A t1       = std::get<2>(lambda_s_t_eff);
  
  A xa(alpha.size(),1.0), tnow(nphon*t1.size(),0.0), tlast(nphon*t1.size(),0.0);
  std::copy( t1.begin(), t1.begin() + t1.size(), tlast.begin() );
  std::copy( t1.begin(), t1.begin() + t1.size(), tnow.begin() );

  F add, exx;
  
  size_t npn = t1.size(), npl = t1.size();

  delta /= tev;

  for( int n = 0; n < nphon; ++n ){
    if ( n > 0 ){ tnow = convol(t1, tlast, delta, npl, npn); }

    for( size_t a = 0; a < alpha.size(); ++a ){
      xa[a] *= lambda_s * alpha[a] * scaling / ( n + 1 );
      exx    = exp(-lambda_s * alpha[a] * scaling)*xa[a];
      for( size_t b = 0; b < beta.size(); ++b ){
        add = exx * interpolate(tnow, beta[b] * sc, betaGrid);
        symSab[b+a*beta.size()] += add < 1e-30 ? 0 : add;
      } // for b in beta
    } // for a in alpha

    npl = npn;

    npn += t1.size() - 1;
    if ( n == 0 ){ continue; }
    if ( npn >= tlast.size() ){ break; }
    for( size_t i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }

  } 
  return std::make_tuple(lambda_s,t_eff);
}



