#include <iostream>
#include <memory>
#include "contin_util/start.h"
#include "contin_util/convol.h"
#include "contin_util/interpolate.h"


auto contin(std::vector<std::vector<std::vector<double>>>& symSab, 
  const std::vector<double>& alpha, const std::vector<double>& beta,
  std::vector<double>& t1, double& delta, const double& tbeta,
  const double& scaling, const double& tev, const double& sc, const int nphon,
  const int itemp ){
  // Start calculates the T1 term, described in Eq. 525 calling start will 
  // also change delta --> delta / tev where tev is temperature in eV. 
  // leapr.f90 calls this deltab
    
  auto lambda_s_t_eff = start( t1, delta, tev, tbeta );
  double lambda_s = std::get<0>(lambda_s_t_eff);

  std::vector<double> xa(alpha.size(),0.0), tnow(nphon*t1.size(),0.0), 
    tlast(nphon*t1.size(),0.0);

  int npn = t1.size();
  std::copy( t1.begin(), t1.begin() + npn, tlast.begin() );

  // This will populate tlast and tnow with blocks, each the same size as t1,
  // corresponding to each order of phonon expansion (nphon). npn is used to 
  // track which block we're at. So we start out with the size of t1, then it
  // will basically go to 2*t1.size(), then 3*t1.size(), etc.

  
  double add, exx;

  // Do the phonon expansion sum 
  // For this, we treat the first iteration slightly different than the all
  // subsequent iterations, because all subsequent iterations require
  // convolution with the one before it. This is following Eq. 526
  
  for( int n = 0; n < nphon; ++n ){

    // If not the first iteration, we'll need to convolve the current t vector
    // with the last one (Following Eq. 526)
    if ( n > 0 ){ tnow = convol(t1, tlast, delta); }

    for( int a = 0; a < alpha.size(); ++a ){
      xa[a] +=  log(lambda_s * alpha[a] * scaling / ( n + 1 ) );

      exx = -lambda_s * alpha[a] * scaling + xa[a];
      // If the exponential value is really small, we'll cut it off since it
      // won't make too much a difference anyway
      if ( exx <= -250.0 ){ continue; }
      exx = exp(exx);

      for( int b = 0; b < beta.size(); ++b ){
        add = n == 0 ? exx * interpolate( t1,   delta, beta[b] * sc ):
                       exx * interpolate( tnow, delta, beta[b] * sc );
        symSab[a][b][itemp] += add < 1e-30 ? 0 : add;
      } // for b in beta
    } // for a in alpha

    if ( n >= 1 ){
      // tnow and tlast will be populated with nphon-many iterations of t1 info,
      // so npn here is being pushed forward by t1 length so that we can get
      // to the next block of vector
      npn += t1.size() - 1;
      for( int i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }
    }
  } // for n in nphon (maxn in leapr.f90) 


  return lambda_s_t_eff;

}

