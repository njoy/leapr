#include <iostream> 
#include <cmath>
#include <vector>
#include "trans_util/s_table_generation.h"
#include "trans_util/sbfill.h"
#include "trans_util/terps.h"


auto trans( const std::vector<double>& alpha, const std::vector<double>& beta,
  const double& trans_weight, double delta, 
  const double& diffusion, const double& sc, const double& arat,
  std::vector<std::vector<std::vector<double>>>& sym_sab,
  const int& itemp, const double& lambda_s ){

  int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;
  std::vector<double> sd(ndmax,0.0);
  std::vector<double> ap(ndmax,0.0);
  std::vector<double> sb(ndmax, 0.0);
  std::vector<double> betan(beta.size(),0.0);

  double nsd;

  // loop over alpha values
  for ( auto a = 0; a < alpha.size(); ++a ){
    double alpha_sc = alpha[a] * sc / arat;
    double ded = 0.4*trans_weight*diffusion*alpha_sc / 
                 sqrt( 1.0 + 1.42*trans_weight*diffusion*diffusion*alpha_sc );
        
    if ( ded == 0 ){ ded = 0.2 * sqrt( trans_weight * alpha_sc );}
    double deb = 10.0 * alpha_sc * delta;
    delta = deb;
    if ( ded < delta ){ delta = ded; }
    nsd = diffusion == 0 ? free_gas_s_table( trans_weight, alpha_sc, ndmax, 
                                             delta, sd, ap ) : 
                           diffusion_s_table( trans_weight, alpha_sc, ndmax, 
                                             delta, sd, ap, diffusion );
    if ( nsd > 1 ){
      for ( int i = 0; i < beta.size(); ++i ){
        betan[i] = beta[i] * sc;
        ap[i] = sym_sab[a][i][itemp];
      }

      // loop over beta values
      for ( int b = 0; b < beta.size(); ++b ){
        double be = betan[b];

        // prepare table of continuous ss on new interval
        b += 1;
        sbfill( sb, nsd, delta, be, ap, betan, b, ndmax );
        b -= 1;

        // convolve s-transport with s-continuous
        double s = 0;
        for ( int i = 0; i < nsd; ++i ){
          double f = 2*(i%2)+2;
          if ( i == 0 or i == nsd - 1 ){ f = 1; }
                    
          s = s + f*sd[i]*sb[nsd+i-1];
          s = s + f*sd[i]*sb[nsd-i-1]*exp(-i*delta);

        }
        s = s < 1e-30 ? 0 : s * delta / 3;
                
        double st = terps(sd,nsd,delta,be);

        if ( st > 0.0 ){ s += exp( -alpha_sc * lambda_s )*st; }
        sym_sab[a][b][itemp] = s;

      } // for beta
    } // if nsd > 0
  } // for alpha
}



