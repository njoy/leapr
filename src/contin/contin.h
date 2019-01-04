#include "contin_util/start.h"
#include "contin_util/convol.h"
#include "contin_util/interpolate.h"
#include "contin_util/checkMoments.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <range/v3/all.hpp>

template <typename A, typename F>
auto contin( const unsigned int itemp, int nphon, F& delta, const F& tbeta, 
  const F& scaling, const F& tev, const F& sc, A t1, const A& alpha, 
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
   * t1      : phonon distribution, (Card12)
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
    
  auto lambda_s_t_eff = start( t1, tbeta, betaGrid );
  F lambda_s = std::get<0>(lambda_s_t_eff),
    t_eff    = std::get<1>(lambda_s_t_eff);

  A xa(alpha.size(),0.0), tnow(nphon*t1.size(),0.0), tlast(nphon*t1.size(),0.0);

  size_t npn = t1.size();
  const size_t np = t1.size();
  std::copy( t1.begin(), t1.begin() + np, tlast.begin() );
  std::copy( t1.begin(), t1.begin() + np, tnow.begin() );

  // This will populate tlast and tnow with blocks, each the same size as t1,
  // corresponding to each order of phonon expansion (nphon). npn is used to 
  // track which block we're at. So we start out with the size of t1, then it
  // will basically go to 2*t1.size(), then 3*t1.size(), etc.

  F add, exx;

  // To be used when checking moments of S(a,b)
  std::vector<int> maxt (1000, 0);
  for ( size_t b = 0; b < beta.size(); ++b ){ maxt[b] = alpha.size() + 1; }

  std::vector<double> eq16(beta.size()), eq14(beta.size());

  // Do the phonon expansion sum 
  // For this, we treat the first iteration slightly different than the all
  // subsequent iterations, because all subsequent iterations require
  // convolution with the one before it. This is following Eq. 526
  
  int npl = np;
  double eq14Val;
  delta /= tev;

  for( int n = 0; n < nphon; ++n ){
    //std::cout << (t1 | ranges::view::all ) << std::endl;
    //std::cout << (tlast | ranges::view::all ) << std::endl;
    //std::cout << std::endl;


    if ( n > 0 ){ tnow = convol(t1, tlast, delta, npl, np, npn); }

    // 
    //if ( n > 0 ){
    //F arat = sc/scaling;
    //auto lambda_s_t_eff_2 = std::make_tuple(lambda_s,t_eff,eq16);
    //return lambda_s_t_eff_2;
    //}
    // 


   
    eq14Val = 0.0;
    for( int b = 0; b < int(beta.size()); ++b ){
      for( int a = 0; a < int(alpha.size()); ++a ){
        if ( b == 0 ) xa[a] +=  log(lambda_s * alpha[a] * scaling / ( n + 1 ) );

        exx = -lambda_s * alpha[a] * scaling + xa[a];
        // make sure to put this in an array or something so that you don't 
        // keep on computing sorry i'm lazy

        if ( exx <= -250.0 ){ continue; }
        exx = exp(exx);
        
        add = exx * interpolate( tnow, delta, beta[b] * sc );
        symSab(a,b,itemp) += add < 1e-30 ? 0 : add;

        if ( symSab(a,b,itemp) != 0        and n >= nphon-1 and 
             add > symSab(a,b,itemp)*0.001 and a < maxt[b] ){ maxt[b] = a; }

      } // for a in alpha
    } // for b in beta

    npl = npn;
    if ( n > 0 ){
      // tnow and tlast will be populated with nphon-many iterations of t1 info,
      // so npn here is being pushed forward by t1 length so that we can get
      // to the next block of vector
      npn += t1.size() - 1;
      if ( npn >= tlast.size() ){ 
        for ( size_t b = 0; b < beta.size(); ++b ){
          double g_prime = 0;
          for ( size_t a = 0; a < alpha.size(); ++a ){ 
            g_prime += symSab(a,b,0); 
          }
          eq14[b] = g_prime; 
          eq16[b] = (b == 0) ? eq14[b] : eq14[b] + eq16[b-1];
        }

      } 

      if ( npn >= tlast.size() ){ 
        auto lambda_s_t_eff_2 = std::make_tuple(lambda_s,t_eff,eq16);
        return lambda_s_t_eff_2; 
      }
      for( size_t i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }
    }
    else { 
      npn += t1.size() - 1;
    }
  } // for n in nphon (maxn in leapr.f90) 

  F arat = sc/scaling;
  //checkMoments( sc, alpha, beta, maxt, itemp, lambda_s, tbeta, arat, t_eff, symSab );

  auto lambda_s_t_eff_2 = std::make_tuple(lambda_s,t_eff,eq16);

  return lambda_s_t_eff_2;

}


