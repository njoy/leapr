#include "translational_util/s_table_generation.h"
#include "translational_util/sbfill.h"
#include "range/v3/all.hpp"
#include "generalTools/tools.h"

template <typename Range, typename Float>
auto translational( Range alpha, Range beta, const Float& transWeight, Float deltaBeta, 
  const Float& diffusion, const Float& lambda_s, const Float& continWgt, 
  Float& t_eff, const Float& temp, Range& sab ){

  using std::exp, std::min; 

  int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;
  Range sabTrans(ndmax), ap(ndmax), sb(ndmax,0.0);

  Float nsd, ded, st, delta, s, f;
  for ( size_t a = 0; a < alpha.size(); ++a ){
    ded = diffusion == 0 ? 
      0.2 * sqrt( transWeight * alpha[a] ) :
      0.4 * transWeight * diffusion * alpha[a] / 
        sqrt( 1.0 + 1.42*transWeight*diffusion*diffusion*alpha[a] );

    delta = min( ded, 10.0 * alpha[a] * deltaBeta );
    nsd = diffusion == 0 ? 
      getFreeGas  ( transWeight, alpha[a], ndmax, delta, sabTrans ) : 
      getDiffusion( transWeight, alpha[a], ndmax, delta, sabTrans, diffusion );

    if (nsd <= 1){ continue; }

    for ( size_t b = 0; b < beta.size(); ++b ){
      ap[b] = sab[beta.size()*a + b];
    }
    for ( size_t b = 0; b < beta.size(); ++b ){
      sbfill( sb, nsd, delta, beta[b], ap, beta, ndmax );
      s = 0;
      for ( int i = 0; i < nsd; ++i ){

        f = (i == 0 or i == nsd - 1) ? 1 : 2*(i%2)+2;
                  
        s += f * sabTrans[i] * sb[nsd+i-1] + 
             f * sabTrans[i] * sb[nsd-i-1] * exp(-i*delta);
      }
      s = (s < 1e-30) ? 0 : s*delta*0.33333333;

      if (beta[b] > nsd*delta){ 
        sab[b+a*beta.size()] = s; 
        continue;
      }
      auto xVals = ranges::view::iota(0,int(sabTrans.size()))
                 | ranges::view::transform([delta](auto x){return delta*x;});
      st = interpolateLog(ranges::view::zip(xVals,sabTrans),beta[b]);

      if ( st > 0.0 ){ s += exp(-alpha[a]*lambda_s)*st; }
      sab[beta.size()*a + b] = s;

    } // for beta
  } // for alpha
  
  t_eff = (continWgt*t_eff + transWeight*temp) / ( continWgt + transWeight );

}

