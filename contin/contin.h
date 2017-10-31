#include <iostream>
#include <memory>
#include "contin_util/start.h"
#include "contin_util/convol.h"
#include "contin_util/interpolate.h"

auto phonon_exp( const double& lambda_s, const double& sc, const double& arat,
                 const double& delta,    const int& nphon, const int& ntempr,
                 const std::vector<double>& alpha, const std::vector<double>& beta,
                 std::vector<double>& t1 ){
 
    std::vector<double> xa(alpha.size(),0.0);
    std::vector<std::vector<std::vector<double>>> sym_sab( 2, 
      std::vector<std::vector<double>> (3, std::vector<double> ( 4, 0.0 ) ) );
    for ( auto entry : sym_sab ){ 
        std::cout << entry.size() << std::endl;
    }
    std::cout << "  " << std::endl;
    std::vector<double> sym_sabOld(alpha.size()*beta.size()*ntempr,0.0);

    double exp_lim = -250.0, tiny = 1e-30;
 
    int npl = t1.size();
    int npn;

    std::vector<int> maxt( 1000, 0 );
    std::vector<double> tnow( nphon*t1.size(), 0.0 );
    std::vector<double> tlast(nphon*t1.size(),0.0);

    for( int i = 0; i < t1.size(); ++i ){ tlast[i] = t1[i]; }

    // Start the phonon expansion with t1
    int c = 0;
    for( int a = 0; a < alpha.size(); ++a ){
        xa[a] = log(alpha[a] * sc * lambda_s / arat );
        double exx = -lambda_s * alpha[a] * sc / arat + xa[a];
        exx = exx > exp_lim ? exp(exx) : 0;
        for( int b = 0; b < beta.size(); ++b ){
            double add = exx * interpolate( t1, delta, beta[b] * sc );
            if ( add < tiny ){ add = 0; }
            sym_sabOld[c] = add;
            c += 1;
        }
    }

    for( int i = 0; i < beta.size(); ++i ){ maxt[i] = alpha.size()+ 1; }
    
    // Do the phonon expansion sum 

    for( int n = 1; n < nphon; ++n ){
        int c = 0;
        npn = t1.size() + npl - 1;
        tnow = convol(t1, tlast, delta, npn, npl);

        for( int a = 0; a < alpha.size(); ++a ){
            xa[a] +=  log(lambda_s * alpha[a] * sc / ( arat * ( n + 1 ) ) );
            double exx = -lambda_s * alpha[a] * sc / arat + xa[a];
            exx = exx > exp_lim ? exp(exx) : 0;
            
            for( int b = 0; b < beta.size(); ++b ){
                double add = exx * interpolate(tnow, delta, beta[b]*sc );
                if ( add < tiny){ add = 0; }
                sym_sabOld[c] += add;
                if (sym_sabOld[c] != 0 and n == nphon-1 and 
                    add > sym_sabOld[c]*.001 and b < maxt[a] ){
                        maxt[a] = b;
                }
                c += 1;
            } // for b in nbeta
        } // for a in nalpha
        npl = npn;
        for( int i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }

    } // for n in nphon (maxn in leapr.f90) 

    return sym_sabOld; 
}



std::vector<double> contin(const int& ntempr, const int& nphon, 
        const std::vector<double>& alpha, const std::vector<double>& beta, 
        int lat, double delta, std::vector<double> phonon_dist, 
        const double& tbeta, const double& arat, const double& tev, 
        const double& sc ){
    
    // Start calculates the T1 term, described in Eq. 525
    // calling start will also change delta --> delta / tev
    // where tev is temperature in eV. leapr.f90 calls this 
    // deltab
    double lambda_s = start( phonon_dist, delta, tev, tbeta );
    std::vector<double> t1 = std::move(phonon_dist);
    return phonon_exp(lambda_s, sc, arat, delta, nphon, ntempr, alpha, beta, t1);

}

