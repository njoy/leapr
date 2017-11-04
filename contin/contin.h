#include <iostream>
#include <memory>
#include "contin_util/start.h"
#include "contin_util/convol.h"
#include "contin_util/interpolate.h"


void phonon_exp( std::vector<std::vector<std::vector<double>>>& sym_sab, 
  const std::vector<double>& alpha, const std::vector<double>& beta,
  std::vector<double>& t1, const double& lambda_s, const double& sc,
  const double& arat, const double& delta, const int nphon, const int itemp ){
 
  std::vector<double> xa(alpha.size(),0.0);

  double exp_lim = -250.0, tiny = 1e-30;
 
  int npn = t1.size();

  std::vector<double> tnow( nphon*t1.size(), 0.0 );
  std::vector<double> tlast(nphon*t1.size(),0.0);

  for( int i = 0; i < t1.size(); ++i ){ tlast[i] = t1[i]; }

  double add, exx;
  // Start the phonon expansion with t1
  for( int a = 0; a < alpha.size(); ++a ){
    xa[a] = log(alpha[a] * sc * lambda_s / arat );
    exx = -lambda_s * alpha[a] * sc / arat + xa[a];
    exx = exx > exp_lim ? exp(exx) : 0;

    for( int b = 0; b < beta.size(); ++b ){
      add = exx * interpolate( t1, delta, beta[b] * sc );
      sym_sab[a][b][itemp] = add < tiny ? 0 : add;

    } // for b in beta
  } // for a in alpha

  // Do the phonon expansion sum 
  for( int n = 1; n < nphon; ++n ){

    npn = t1.size() + npn - 1;
    tnow = convol(t1, tlast, delta);

    for( int a = 0; a < alpha.size(); ++a ){
      xa[a] +=  log(lambda_s * alpha[a] * sc / ( arat * ( n + 1 ) ) );
      exx = -lambda_s * alpha[a] * sc / arat + xa[a];
      exx = exx > exp_lim ? exp(exx) : 0;
            
      for( int b = 0; b < beta.size(); ++b ){
        add = exx * interpolate(tnow, delta, beta[b]*sc );
        sym_sab[a][b][itemp] += add < tiny ? 0 : add;

      } // for b in beta
    } // for a in alpha

    for( int i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }
  
  } // for n in nphon (maxn in leapr.f90) 
}



auto contin(std::vector<std::vector<std::vector<double>>>& sym_sab, 
  const std::vector<double>& alpha, const std::vector<double>& beta,
  std::vector<double>& phonon_dist, double& delta, const double& tbeta,
  const double& arat, const double& tev, const double& sc, const int nphon,
  const int itemp ){

  // Start calculates the T1 term, described in Eq. 525 calling start will 
  // also change delta --> delta / tev where tev is temperature in eV. 
  // leapr.f90 calls this deltab
    
  double lambda_s = start( phonon_dist, delta, tev, tbeta );
  phonon_exp( sym_sab, alpha, beta, phonon_dist, lambda_s, sc, arat, delta, 
              nphon, itemp );
  return lambda_s;

}

