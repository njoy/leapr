
#ifndef THERMR_SIG
#define THERMR_SIG


#include "calcem/calcem_util/sig_util/terpq.h"
#include "general_util/sigfig.h"
#include <cmath>

inline double do160(double sigc, double sb, double s, double bb, double sabflg,
  double sigmin ){
  // Based off of Eq. 225, it appears that s = S(a,b). 
  // sigc = sqrt(E'/E)/(2 kb T), sb = sigma_b, and bb = beta. 
  // Note that calling do160 at the very and of sig seems mysterious, because
  // S(a,b) = S(a0,b0) + 1/2 ln(a0/a) - a0/a * b^2/b0^2 ( S(a0,b0) - S(a0,b1) )
  // and I really don't know where this comes from. (note the last term is b/c
  // of variable cliq that is calculated elsewhere and fed into sig)
  if (s-bb/2 <= sabflg ){ return 0.0; }
  double sigVal = sigc * sb * exp(s-bb/2);
  return (sigVal < sigmin) ? 0 : sigVal;
}

inline auto do155( int ia, int ib, const std::vector<double>& alpha,
    const std::vector<double>& beta, const std::vector<std::vector<double>>& sab,
    double a, double bbb, double sigc, double sb, double bb, double sabflg,
    double sigmin ){
   if ( (unsigned) ia+1 == alpha.size()) ia -= 1;
   if ( (unsigned) ib+1 == beta.size())  ib -= 1;
   ia -= 1;
   ib -= 1;
   double s1 = terpq( alpha[ia],     alpha[ia+1],     alpha[ia+2], a, 
                      sab[ia][ib],   sab[ia+1][ib],   sab[ia+2][ib] );
   double s2 = terpq( alpha[ia],     alpha[ia+1],     alpha[ia+2], a, 
                      sab[ia][ib+1], sab[ia+1][ib+1], sab[ia+2][ib+1] );
   double s3 = terpq( alpha[ia],     alpha[ia+1],     alpha[ia+2], a, 
                      sab[ia][ib+2], sab[ia+1][ib+2], sab[ia+2][ib+2] );
   double s  = terpq( beta[ib], beta[ib+1], beta[ib+2], bbb, s1, s2, s3 );
   return do160(sigc, sb, s, bb, sabflg, sigmin );
}


inline auto doSCTApproximation( double a, double teff, 
  double sigc, 
  double sb2, double tev, double sigmin, double sabflg, double bb, 
  double s, double sb ){
 /* The SCT approximation is calculated, according to Eq. 230. Note that since
  * some evaluations give S(a,b) for a molecule or compund (e.g. C6H6 or BeO), 
  * the correspondint SCT approximation must contin terms for both atoms.
  * So in this case (when s2 > 0) we need to recalculate Eq. 230 taking into 
  * account the secondary atom's parameters.
  */

  // short collision time for large beta or alpha

  // Eq. 230 in the NJOY manual
  double arg, sigVal = 0;
  for ( const double& sb_val : {sb, sb2} ){
    // Eq. 230 in the NJOY manual
    arg = (a-std::abs(bb))*(a-std::abs(bb))*tev/(4.0*a*teff) + (std::abs(bb)+bb)/2.0;
    s = (-arg > sabflg) ? exp(-arg)/(sqrt(4.0*M_PI*a*teff/tev)) : 0;
    sigVal += sigc * sb_val * s;
  }
  // If sb2 > 0, then we're going to have to add the second atom's contribution
  // to sig. 

  return sigVal < sigmin ? 0 : sigVal;

}


inline auto sig( const double& e, const double& ep, const double& u, 
  const double& tev,  
  const std::vector<double>& alpha, const std::vector<double>& beta,
  const std::vector<std::vector<double>>& sab, const double az, 
  const double tevz, const int lasym,// const double az2, const double teff2, 
  const int lat, const double cliq, const double sb, const double sb2, 
  const double teff, const int iinc ){

 /*-------------------------------------------------------------------
  * Compute the differential scattering cross section from e to
  * ep through the angle with the cosine u from endf tabulated
  * data or an analytic law.
  *-------------------------------------------------------------------
  */
  int i,ib,ia;
  double bb,a,sigc,b,s,s1,s2,s3,bbm;
  double sigmin=1.e-10, sabflg=-225.e0, amin=1.e-6, test1=0.2e0,
         test2=30.e0;
  double sigVal;

  // common factors.
  sigc=sqrt(ep/e)/(2.0*tev);
  double a_tev   = (e+ep-2*u*sqrt(e*ep))/(az*tev);
  double a_tevz  = (e+ep-2*u*sqrt(e*ep))/(az*tevz);
  double bb_tev  = (ep-e)/tev;     
  double bb_tevz = (ep-e)/tevz;
  if ( a_tev < amin ){ a_tev = amin; a_tevz = amin * tev / tevz; }


  // free-gas scattering. Plug Eq. 229 into Eq. 225, solve.
  if ( iinc == 1 ){ 

    // This computes S(a,b) according to Eq. 229.
    // S(a,b) = 1/sqrt(4pi*a) * exp( - ( a^2 + b^2 ) / 4a )          [Eq. 229]
    double sab = 1.0/sqrt(4.0*M_PI*a_tev) * exp( -(a_tev*a_tev + bb_tev*bb_tev)/(4.0*a_tev) );
    // sigma = sigma_b / 2 kb T  * sqrt(E'/E) * exp( -b/2 ) * S(a,b) [Eq. 225]
    // Note that sigma_b is the bound scattering cross section, which is 
    // defined by Eq. 228 to be 
    //             sigma_b = sigma_f * ( A + 1 )^2 / A^2             [Eq. 228]
    // where sigma_f is the characteristic free scattering cross section
    sigVal = sb/(2.0*tev) * sqrt(ep/e) * exp(-bb_tev/2.0) * sab;
    return sigVal < sigmin ? 0 : sigVal;

  } // end free gas scattering option

  if (iinc != 2){
     // other options not yet implemented.
     throw std::exception(); // call error('sig','illegal option.',' ')
   }

  // Now we change the definitions of a and b so that they're defined using
  // tevz instead of tev
  if ( lat == 1 ){ 
    b = std::abs(bb_tevz); a = a_tevz;
  }
  else {
    b = std::abs(bb_tev); a = a_tev;
  }
  //std::cout << a << std::endl;
  b = sigfig(b,8,0);
  //if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "in sig     " << b <<  std::endl;

  if (a > alpha[alpha.size()-1]) {
    // go to 170
    return doSCTApproximation( a_tev, teff, 
        sigc, sb2, tev, sigmin, 
        sabflg, bb_tev, s, sb );
  }

  if (lasym == 1) {
     // S(a,b) is not symmetric, probably because we have orthohydrogen or 
     // parahydrogen.
     bbm = (lat == 1) ? bb_tevz : bb_tev;
     if ( bbm > beta[beta.size()-1] or bbm < beta[0] ){
       // go to 170
       return doSCTApproximation( a_tev, teff, 
           sigc, sb2, tev, 
           sigmin, sabflg, bb_tev, s, sb );
     } 
  } // end if 
  else {
    if ( b > beta[beta.size()-1] ){
       return doSCTApproximation( a_tev, teff, 
           sigc, sb2, tev, 
           sigmin, sabflg, bb_tev, s, sb );
     }
  } // end if

  if (cliq == 0.0 or a >= alpha[0] or lasym == 1 or b >= test1 ){ 

    // Find index ia such that alpha[ia-1] < a < alpha[ia]. Same for ib
    auto ia1 = std::distance(alpha.cbegin(), std::lower_bound(alpha.cbegin(), 
          alpha.cend(), a));
    auto ib1 = std::distance(beta.cbegin(),  std::lower_bound(beta.cbegin(), 
          beta.cend(), b));
    ia = ia1 > 1 ? ia1 : 1;
    ib = ib1 > 1 ? ib1 : 1;
    //std::cout << "-----  " << ia << "   " << ib  << std::endl;
    //std::cout << "-----  " << a << "   " << b  << std::endl;

    // 150 continue
    if (a*az < test2 and b < test2) { 
      return do155( ia, ib, alpha, beta, sab, a, b, sigc, sb,
          bb_tev, sabflg, sigmin );
    }
    if ( a*az >= test2 or b >= test2 ){
      if ( sab[ia-1][ib-1]   <= sabflg or sab[ia+1-1][ib-1]   <= sabflg or 
           sab[ia-1][ib+1-1] <= sabflg or sab[ia+1-1][ib+1-1] <= sabflg ) {
        return doSCTApproximation( a_tev, teff, 
         sigc, sb2, tev, 
          sigmin, sabflg, bb_tev, s, sb );
      } 
      else { 
        return do155( ia, ib, alpha, beta, sab, a, b, sigc, sb,
          bb_tev, sabflg, sigmin );
      }
    }
  }


  s=sab[0][0]+log(alpha[0]/a)/2-cliq*b*b/a;
  if (s < sabflg) s=sabflg;
  // go to 160
  return do160(sigc, sb, s, bb_tev, sabflg, sigmin );

}



#endif
