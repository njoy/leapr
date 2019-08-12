#include "trans_util/s_table_generation.h"
#include "trans_util/sbfill.h"
#include "trans_util/terps.h"
#include "range/v3/all.hpp"

template <typename Range, typename Float>
auto trans( Range alpha, Range beta, const Float& transWeight, 
  Float delta, const Float& diffusion, const Float& sc, const Float& scaling, 
  const Float& lambda_s, const Float& tbeta, Float& t_eff, 
  const Float& temp, Range& sym_sab ){
  using std::exp, std::min; 

  for (auto& a : alpha ){ a *= scaling; }
  for (auto& b : beta  ){ b *= sc;      }

  int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;
  Range sabTrans(ndmax), ap(ndmax), sab(ndmax);

  Float nsd, ded, st;
  for ( size_t a = 0; a < alpha.size(); ++a ){
    ded = diffusion == 0 ? 
      0.2 * sqrt( transWeight * alpha[a] ) :
      0.4 * transWeight * diffusion * alpha[a] / 
        sqrt( 1.0 + 1.42*transWeight*diffusion*diffusion*alpha[a] );

    delta = min( ded, 10.0 * alpha[a] * delta );
    nsd = diffusion == 0 ? 
      getFreeGas  ( transWeight, alpha[a], ndmax, delta, sabTrans ) : 
      getDiffusion( transWeight, alpha[a], ndmax, delta, sabTrans, diffusion );

    if ( nsd > 1 ){
      for ( size_t b = 0; b < beta.size(); ++b ){
        ap[b] = sym_sab[beta.size()*a + b];
      }

      for ( size_t b = 0; b < beta.size(); ++b ){
        Float be = beta[b];
        sbfill( sab, nsd, delta, be, ap, beta, ndmax );
        Float s = 0;
        for ( int i = 0; i < nsd; ++i ){
          Float f = 2*(i%2)+2;
          if ( i == 0 or i == nsd - 1 ){ f = 1; }
                    
          s += f * sabTrans[i] * sab[nsd+i-1] + 
               f * sabTrans[i] * sab[nsd-i-1] * exp(-i*delta);

        }
        s = (s < 1e-30) ? 0 : s*delta*0.33333333;

        st = terps(sabTrans,delta,be);
        if ( st > 0.0 ){ s += exp(-alpha[a]*lambda_s)*st; }

        sym_sab[beta.size()*a + b] = s;

      } // for beta
    } // if nsd > 0
  } // for alpha
  
  t_eff = (tbeta*t_eff + transWeight*temp) /
                     ( tbeta + transWeight );

}

