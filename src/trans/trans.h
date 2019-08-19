#include "trans_util/s_table_generation.h"
#include "trans_util/sbfill.h"
#include "trans_util/terps.h"
#include "range/v3/all.hpp"

template <typename Range, typename Float>
auto trans( Range alpha, Range beta, const Float& transWeight, Float deltaBeta, 
  const Float& diffusion, const Float& sc, const Float& scaling, 
  const Float& lambda_s, const Float& tbeta, Float& t_eff, const Float& temp, 
  Range& existingSAB ){

  using std::exp, std::min; 

  for (auto& a : alpha ){ a *= scaling; }
  for (auto& b : beta  ){ b *= sc;      }

  int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;
  Range sabTrans(ndmax), ap(ndmax), sab(ndmax);

  Float nsd, ded, st, delta;
  for ( size_t a = 0; a < alpha.size(); ++a ){
    ded = diffusion == 0 ? 
      0.2 * sqrt( transWeight * alpha[a] ) :
      0.4 * transWeight * diffusion * alpha[a] / 
        sqrt( 1.0 + 1.42*transWeight*diffusion*diffusion*alpha[a] );

    delta = min( ded, 10.0 * alpha[a] * deltaBeta );
    nsd = diffusion == 0 ? 
      getFreeGas  ( transWeight, alpha[a], ndmax, delta, sabTrans ) : 
      getDiffusion( transWeight, alpha[a], ndmax, delta, sabTrans, diffusion );
    //std::cout << a+1 << "    " << "nsd   " << nsd << std::endl;
    //std::cout << a+1 << "    " << sabTrans[0] << "   " << sabTrans[1] << "   " << sabTrans[2] << std::endl;

    if ( nsd > 1 ){
      for ( size_t b = 0; b < beta.size(); ++b ){
        ap[b] = existingSAB[beta.size()*a + b];
      }
    //std::cout << a+1 << "    " << ap[0] << "   " << ap[1] << "   " << ap[2] << std::endl;

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
        //std::cout << a+1 << "    " << b+1 << "     " << s << std::endl;
        s = (s < 1e-30) ? 0 : s*delta*0.33333333;

        st = terps(sabTrans,delta,be,nsd);
        if ( st > 0.0 ){ s += exp(-alpha[a]*lambda_s)*st; }

        existingSAB[beta.size()*a + b] = s;

      } // for beta
    }
  } // for alpha
  
  t_eff = (tbeta*t_eff + transWeight*temp) /
                     ( tbeta + transWeight );

}

