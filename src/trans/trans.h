#include "trans_util/s_table_generation.h"
#include "trans_util/sbfill.h"
#include "trans_util/terps.h"
#include "range/v3/all.hpp"


template <typename Range, typename Float, typename Range2>
auto trans( const Range& alpha, const Range& beta, const Float& trans_weight, 
  Float delta, const Float& diffusion, const Float& sc, const Float& scaling, 
  const int& itemp, const Float& lambda_s, const Float& tbeta, Range& t_eff_vec, 
  const Range& temp_vec, Range2& sym_sab ){
  using std::exp; 

  Float deltaInitial = delta;
  int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;
  Range sabTrans(ndmax), ap(ndmax), sab(ndmax), betan(beta.size());

  Float nsd, alpha_sc, ded;
  for ( size_t a = 0; a < alpha.size(); ++a ){
    alpha_sc = alpha[a] * scaling;

    ded = diffusion == 0 ? 
      0.2 * sqrt( trans_weight * alpha_sc ) :
      0.4 * trans_weight * diffusion * alpha_sc / 
        sqrt( 1.0 + 1.42*trans_weight*diffusion*diffusion*alpha_sc );

    delta = std::min( ded, 10.0 * alpha_sc * deltaInitial );
    nsd = diffusion == 0 ? getFreeGas ( trans_weight, alpha_sc, ndmax, 
                                             delta, sabTrans ) : 
                           getDiffusion( trans_weight, alpha_sc, ndmax, 
                                             delta, sabTrans, diffusion );
    if ( nsd > 1 ){
      for ( size_t b = 0; b < beta.size(); ++b ){
        betan[b] = beta[b] * sc;
        ap[b] = sym_sab[beta.size()*a + b];
      }

      for ( size_t b = 0; b < beta.size(); ++b ){
        Float be = betan[b];
        sbfill( sab, nsd, delta, be, ap, betan, ndmax );
        Float s = 0;
        for ( int i = 0; i < nsd; ++i ){
          Float f = 2*(i%2)+2;
          if ( i == 0 or i == nsd - 1 ){ f = 1; }
                    
          s += f * sabTrans[i] * sab[nsd+i-1] + 
               f * sabTrans[i] * sab[nsd-i-1] * exp(-i*delta);

        }
        s = (s < 1e-30) ? 0 : s*delta*0.33333333;

        Float st = terps(sabTrans,delta,be);

        if ( st > 0.0 ){ s += exp(-alpha_sc*lambda_s)*st; }
        sym_sab[beta.size()*a + b] = s;

      } // for beta
    } // if nsd > 0
  } // for alpha
  
  t_eff_vec[itemp] = (tbeta*t_eff_vec[itemp] + trans_weight*temp_vec[itemp]) /
                     ( tbeta + trans_weight );

}

