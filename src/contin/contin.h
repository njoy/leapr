#include "contin_util/start.h"
#include "contin_util/convol.h"
#include "contin_util/interpolate.h"
#include "contin_util/checkMoments.h"
#include <unsupported/Eigen/CXX11/Tensor>

template <typename A, typename F>
auto contin( const unsigned int itemp, int nphon, F& delta, const F& tbeta, 
  const F& scaling, const F& tev, const F& sc, A rho, const A& alpha, 
  const A& beta, Eigen::Tensor<F,3>& symSab, A betaGrid ){

  /* Inputs
   * ------------------------------------------------------------------------
   * itemp   : temperature index, to help with updating the S(a,b) vector with
   *           final result
   * nphon   : phonon expansion order, an input provided by Card3. This decides 
   *           when to stop the sum in Eq. 522. 
   * tbeta   : normalization for the continuous part. This decides how heavily
   *           weighted the solid type spectrum will be weighted. Card13
   * scaling : sc / arat, where sc is a user-influenced scaling term and arat
   *           is a mass ratio. Computed in leapr.cpp and applied to all a,b.
   * tev     : temperature in eV. Computed in leapr.cpp vie T * k_b.
   * sc      : user influenced scaling term
   * rho     : phonon distribution, (Card12)
   * alpha   : alpha values (Card 8)
   * beta    : beta values  (Card 9)
   * symSab  : symmetric S(a,b). S[1][2][3] = S(a=1,b=2,T=3). Starts blank.
   *    
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Calls start to calculate the Debye-Waller coefficient (lambda), the 
   *   effective temperature, and T_1. These values are defined by Eq. 521,
   *   530, and 525 respectively.
   * * Start solving Eq. 522 - 523. The xa vector accounts for the 
   *            exp(-lambda*alpha) * ( alpha * lambda )^n / n!
   *   contribution. tnow is T_n, tlast is T_n-1. Convol is used to advance 
   *   from a given T_n to the next T_n+1, and follows Eq. 526.
   * * After Eq. 522 - 523 is solved, record it in S(a,b) vector. 
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * SymSab is the modified S(a,b), following Eq. 522 - 523
   */


  // Start calculates the T1 term, described in Eq. 525 calling start will 
  // also change delta --> delta / tev where tev is temperature in eV. 
  // leapr.f90 calls this deltab
  
  for ( F& x : betaGrid ){ x /= tev; }
    
  auto lambda_s_t_eff = start( rho, tbeta, betaGrid );
  F lambda_s = std::get<0>(lambda_s_t_eff),
    t_eff    = std::get<1>(lambda_s_t_eff);

  auto t1_temp = std::get<2>(lambda_s_t_eff);

  A t1(t1_temp.size());

  size_t npn = t1.size();
  const size_t np = t1.size();
  std::copy( t1_temp.begin(), t1_temp.begin() + np, t1.begin() );

  A xa(alpha.size(),1.0), tnow(nphon*t1.size(),0.0), tlast(nphon*t1.size(),0.0);
  std::copy( t1.begin(), t1.begin() + np, tlast.begin() );
  std::copy( t1.begin(), t1.begin() + np, tnow.begin() );

  // This will populate tlast and tnow with blocks, each the same size as t1,
  // corresponding to each order of phonon expansion (nphon). npn is used to 
  // track which block we're at. So we start out with the size of t1, then it
  // will basically go to 2*t1.size(), then 3*t1.size(), etc.

  F add, exx, g_prime;

  // Do the phonon expansion sum 
  // For this, we treat the first iteration slightly different than the all
  // subsequent iterations, because all subsequent iterations require
  // convolution with the one before it. This is following Eq. 526
  
  int npl = np;
  delta /= tev;

  for( int n = 0; n < nphon; ++n ){
//    if ( n > 0 ){ tnow = convol(t1, tlast, npn, betaGrid); }
    if ( n > 0 ){ tnow = convol(t1, tlast, delta, npl, np, npn); }

    for( size_t a = 0; a < alpha.size(); ++a ){
      xa[a] *=  lambda_s * alpha[a] * scaling / ( n + 1 );
      exx = exp(-lambda_s * alpha[a] * scaling)*xa[a];
      for( size_t b = 0; b < beta.size(); ++b ){
        add = exx * interpolate( tnow, beta[b] * sc,betaGrid );
        symSab(a,b,itemp) += add < 1e-30 ? 0 : add;
      } // for b in beta
    } // for a in alpha

    npl = npn;

    // tnow and tlast will be populated with nphon-many iterations of t1 info,
    // so npn here is being pushed forward by t1 length so that we can get
    // to the next block of vector
    npn += t1.size() - 1;
    if ( n == 0 ){ 
      continue;
    }

    
    if ( npn >= tlast.size() ){ 
      return std::make_tuple(lambda_s,t_eff); }

    for( size_t i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }

  } // for n in nphon (maxn in leapr.f90) 

  //for ( auto x : alpha0Additions ) { std::cout << x << std::endl; }
  return std::make_tuple(lambda_s,t_eff);
}


