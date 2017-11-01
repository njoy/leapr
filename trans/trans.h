#include <iostream> 
#include <cmath>
#include <vector>
#include "trans_util/s_table_generation.h"
#include "trans_util/sbfill.h"
#include "trans_util/terps.h"


auto trans( const std::vector<double>& alpha, const std::vector<double>& beta,
            const int lat, const double& trans_weight, double delta, 
            const double& diffusion_const, const double& sc, const double& arat,
            const double& tev, std::vector<std::vector<std::vector<double>>>& ssm,
            const int& itemp, const double& lambda_s, 
            std::vector<std::vector<std::vector<double>>>& sym_sab ){
    double c0 = 0.4;    double c1 = 1.0;    double c2 = 1.42;
    double c3 = 0.2;    double c4 = 10.0;

    int ndmax = beta.size() > 1e6 ? beta.size() : 1e6;
    std::vector<double> sd(ndmax,0.0);
    std::vector<double> ap(ndmax,0.0);
    std::vector<double> sb(ndmax, 0.0);
    std::vector<double> betan(beta.size(),0.0);
    double nsd = 0;

    for ( auto a = 0; a < alpha.size(); ++a ){
        double alpha_sc = alpha[a] * sc / arat;
        double ded = c0 * trans_weight * diffusion_const * alpha_sc / 
                 sqrt( c1 + c2 * trans_weight * diffusion_const * alpha_sc );

        if ( ded == 0 ){ ded = c3 * sqrt( trans_weight * alpha_sc );}
        double deb = c4 * alpha_sc * delta;
        delta = deb;
        if ( ded < delta ){ delta = ded; }
        if (diffusion_const == 0){
            free_gas_s_table( trans_weight, alpha_sc, ndmax, delta, sd, ap, nsd );
        }
        else{
            diffusion_s_table( diffusion_const, trans_weight, alpha_sc, ndmax, 
            delta, sd, ap, nsd );
        }

        if ( nsd > 1 ){
            for ( int i = 0; i < beta.size(); ++i ){
                betan[i] = beta[i] * sc;
                ap[i] = ssm[a][i][itemp];
            }
            // loop over beta values
            for ( int ibeta = 0; ibeta < beta.size(); ++ibeta ){
                double s = 0;
                double be = betan[ibeta];
                // prepare table of continuous ss on new interval
                int nbt = nsd;
                ibeta += 1;
                sbfill( sb, nbt, delta, be, ap, betan, beta.size(), ibeta, ndmax );
                ibeta -= 1;
                // convolve s-transport with s-continuous
                for ( int i = 0; i < nbt; ++i ){
                    double f = 2*(i%2)+2;
                    if ( i == 0 or i == nbt - 1 ){ f = 1; }
                    s = s + f*sd[i]*sb[nbt+i-1];
                    s = s + f*sd[i]*sb[nbt-i+1]*exp(-i*delta);
                }
                s = s*delta/3;
                if ( s < 1e-30 ){ s = 0; }
                
                double st = terps(sd,nbt,delta,be);
                if ( st > 0.0 ){ s += exp( -alpha_sc * lambda_s )*st; }
                sym_sab[a][ibeta][itemp] = s;

            }
        }
    }
}



