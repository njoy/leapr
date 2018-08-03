#include "contin_util/start.h"
#include "contin_util/convol.h"
#include "contin_util/interpolate.h"
#include "contin_util/checkMoments.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include <range/v3/all.hpp>

template <typename A, typename F>
auto contin( const unsigned int itemp, int nphon, F& delta, const F& tbeta, 
  const F& scaling, const F& tev, const F& sc, A t1, const A& alpha, 
  const A& beta, Eigen::Tensor<F,3>& symSab ){

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
    
  auto lambda_s_t_eff = start( t1, delta, tev, tbeta );
  F lambda_s = std::get<0>(lambda_s_t_eff),
    t_eff    = std::get<1>(lambda_s_t_eff);

  A xa(alpha.size(),0.0), tnow(nphon*t1.size(),0.0), 
    tlast(nphon*t1.size(),0.0);

  size_t npn = t1.size();
  std::copy( t1.begin(), t1.begin() + npn, tlast.begin() );
  std::copy( t1.begin(), t1.begin() + npn, tnow.begin() );

  // This will populate tlast and tnow with blocks, each the same size as t1,
  // corresponding to each order of phonon expansion (nphon). npn is used to 
  // track which block we're at. So we start out with the size of t1, then it
  // will basically go to 2*t1.size(), then 3*t1.size(), etc.

  F add, exx;

  // To be used when checking moments of S(a,b)
  std::vector<int> maxt (1000, 0);
  for ( size_t b = 0; b < beta.size(); ++b ){
    maxt[b] = alpha.size() + 1;
  }

  // Do the phonon expansion sum 
  // For this, we treat the first iteration slightly different than the all
  // subsequent iterations, because all subsequent iterations require
  // convolution with the one before it. This is following Eq. 526
  
  /*
  int counter = 0;
  std::vector<double> correctXa {-3.1687178586617719, -7.0305828978834892, -11.297913045213370, -15.852925264995033, -20.631081036090904 };
  std::vector<double> correctEx {-3.2107753456537758, -7.0726403848754931, -11.339970532205374, -15.894982751987037, -20.673138523082908 };
  std::vector<double> correctExx { 4.0325335074603723E-002, 8.4799112767416690E-004, 1.1888125273830262E-005, 1.2499616851585645E-007, 1.0514049462811860E-009 };
  std::vector<double> correctSsm500 { 4.3093012985286720e-2, 4.3615278228135219e-2, 4.3620924874659871e-2, 4.3620975001566170e-2, 4.3620975370788316e-2 };
  */
  std::vector<double> correctSsm050 { 0.23972713036544355, 0.26219351437564831, 0.26397184441734067, 0.26408986955492991, 0.26409638252188933, 0.26409668821687571 };
  std::vector<double> correctAdd { 0.23972713036544355, 2.2466384010204784E-002, 1.7783300416923651E-003, 1.1802513758923018E-004, 6.5129669594191189E-006, 3.0569498636212804E-007 };
  std::vector<double> correctXa { -1.1538148381195072, -3.0007768567989599, -5.2532039835865767, -7.7933131828259743, -10.556565933379581, -13.502140240727144 };
  std::vector<double> correctExx {0.23009891654938383, 3.6290183211186905E-002, 3.8156847708548046E-003, 3.0089646115464809E-004, 1.8982423501427329E-005, 9.9794128685997906E-007 };
  std::vector<double> correctSt { 1.0418438033539950, 0.61907607022714728, 0.46605790270614433, 0.39224501722727217, 0.34310513401670734, 0.30632562294721455 };

  for( int n = 0; n < nphon; ++n ){
    if ( n > 0 ){ tnow = convol(t1, tlast, delta); }
   
    for( int a = 0; a < int(alpha.size()); ++a ){
      xa[a] +=  log(lambda_s * alpha[a] * scaling / ( n + 1 ) );

      exx = -lambda_s * alpha[a] * scaling + xa[a];

      //if ( exx <= -250.0 ){ continue; }
      if ( exx <= -250.0 ){ exx = 0.0; }
      else { exx = exp(exx); }
        
      for( int b = 0; b < int(beta.size()); ++b ){
        add = exx * interpolate( tnow, delta, beta[b] * sc );
        symSab(a,b,itemp) += add < 1e-30 ? 0 : add;

        if ( symSab(a,b,itemp) != 0 and n >= nphon-1 ) {
          if (add > symSab(a,b,itemp)*0.001 and a < maxt[b]){ maxt[b] = a; }
        } 
        //if ( b == 5 and a == 0 ){ std::cout << std::setprecision(15) << n+1 << "     " << add << "     " << std::abs(symSab(a,b,0) - correctSsm500[counter++]) << std::endl; }
        if ( a == 5 and b == 0 and n == 3 ){ std::cout << std::setprecision(15) << n+1 << "     " << tnow[0] << "    " << tnow[1] << "    " << tnow[2] << "     " << tnow[3]<< std::endl; }
        //if ( a == 5 and b == 0 ){ std::cout << std::setprecision(15) << n+1 << "     " << add << "     " << symSab(a,b,0) << "    " << correctSsm050[counter++]<< std::endl; }

      } // for b in beta
      //if ( a == 0 ){ std::cout << std::setprecision(16) << n+1 << "    " << std::abs(symSab(0,5,0)-correctSsm500[counter++]) << std::endl; }
    } // for a in alpha

    if ( n > 0 ){
      // tnow and tlast will be populated with nphon-many iterations of t1 info,
      // so npn here is being pushed forward by t1 length so that we can get
      // to the next block of vector
      npn += t1.size() - 1;
      for( size_t i = 0; i < npn; ++i ){ tlast[i] = tnow[i]; }
    }
    else { 
      npn += t1.size() - 1;
    }
  } // for n in nphon (maxn in leapr.f90) 

  F arat = sc/scaling;
  checkMoments( sc, alpha, beta, maxt, itemp, lambda_s, tbeta, arat, t_eff, symSab );

  return lambda_s_t_eff;

}

