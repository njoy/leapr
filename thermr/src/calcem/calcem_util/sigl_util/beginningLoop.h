#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"

double maxOf3Vals( const double& a, const double& b, const double& c ){
  return (a < b) ? (b < c ? c : b) : (a < c ? c : a);
}

auto adaptiveLinearization( std::vector<double>& x, std::vector<double>& y, 
  const double& e, const double& ep, const double& tev, const double& tevz, 
  const std::vector<double>& alpha, const std::vector<double>& beta,
  const std::vector<std::vector<double>>& sab, 
  const double& az, const double& az2, const int& lasym, const double& teff, 
  const double& teff2, const int& lat, const double& cliq, const double& sb, 
  const double& sb2, const int& iinc, 
  double& xl, const double& eps, const double& seep, 
  const double& s1bb  ){
  /* So here we consider three angles - a cosine value of -1, a cosine value 
   * that corresponds to alpha = sqrt(1+beta^2) [According to Eq. 227], and a
   * cosine value of 1. We calculate the incoherent cross sections for all 
   * of these angles (using sig), and determine the maximum cross section of 
   * these three.
   */

  // prime stack for equally-probable angles
  //std::cout << 130 << std::endl;
  // adaptive linearization
  // Consider a cosine mu equal to -1. What's the cross section?
  x[2] = -1;
//  if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "going into sig     " << y[2] << std::endl;
  y[2] = sig(e,ep,x[2],tev,alpha,beta,sab,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
//  if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "finished sig     " << y[2] << std::endl;
//  if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 ) return 0.0;

  // Consider a cosine mu that corresponds to an alpha value of sqrt(1+beta^2).
  // What's the cross section?
  x[1] = 0.5 * seep * ( e + ep - (s1bb-1) * az * tev );
  if (abs(x[1]) > 1-eps) x[1] = 0.99;
  x[1] = sigfig(x[1],8,0);
  y[1] = sig(e,ep,x[1],tev,alpha,beta,sab,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  // Consider a cosine mu equal to 1. What's the cross section?
  x[0] = 1;
  y[0] = sig(e,ep,x[0],tev,alpha,beta,sab,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  double ymax = maxOf3Vals(y[0],y[1],y[2]);
  //if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "in sigl     " << y[1] << "     " << y[2] << std::endl;
  return ( ymax < eps ) ? eps : ymax;

}



void shiftOver( int& i, std::vector<double>& x, std::vector<double>& y, 
  double& xm, double& yt ){
  // x = [   mu1      mu2          mu3            0    0   0 ... ]
  // y = [ s(mu1)   s(mu2)       s(mu3)           0    0   0 ... ]
  //   |
  //   |  
  //   V
  // x = [   mu1      mu2      0.5*(mu2+mu3)     mu3   0   0 ... ]
  // y = [ s(mu1)   s(mu2)   s(0.5*(mu2+mu3))  s(mu3)  0   0 ... ]
 
  i = i + 1;
  x[i-1] = x[i-2];
  y[i-1] = y[i-2];
  x[i-2] = xm;
  y[i-2] = yt;
}

auto do_110(int& i, std::vector<double>& x, std::vector<double>& y, 
  const double& e, const double& ep, const double& tev, 
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<std::vector<double>>& sab, const 
  double& az, const double tevz, const int& lasym, const double& az2, const 
  double& teff2, const int& lat, const double& cliq, const double& sb, const 
  double& sb2, const double& teff, const int& iinc, const double& xtol, 
  const double& tol, double& ymax){
  /* So note that here, x[i-2] = mu_a, x[i-1] = mu_c. Similarly, y[i-2] is 
   * the incoh. cross section [Eq. 225] corresponding to mu_a, and y[i-1] is
   * the incoh. cross section [Eq. 225] corresponding to mu_c.
   *
   * s(mu_1) is the incoh. cross section corresponding to cosine mu_1
   *
   * We check to see if the grid is fine enough. 
   * 
   * Consider mu_b = 1/2 * ( mu_a + mu_c )
   * If we trust our grid, then 
   *  
   *   s(mu_b) should be adequately described by 1/2 * ( s(mu_a) + s(mu_c) )
   *
   * So we compute this approximation (ym), and calculate the real deal. 
   * If they're not close enought, then bisect the grid spacing and spread
   * things out more (visualization below).
   *
   */

  double xm, ym, yt;
  while ( i < x.size() ){ 
    // std::cout << 110 << std::endl;
    
    xm = 0.5*( x[i-2] + x[i-1] );
    xm = sigfig(xm,8,0);
    ym = 0.5*( y[i-2] + y[i-1] );
    yt = sig(e, ep, xm, tev, alpha, beta, sab, az, tevz, lasym, az2, teff2,
      lat, cliq, sb, sb2, teff, iinc);
    
    if ( ( abs(yt-ym) <= tol*abs(yt)+tol*ymax/50.0 and 
           abs(y[i-2]-y[i-1]) <= ym+ymax/100.0 and 
           (x[i-2]-x[i-1]) < 0.5 ) or
         ( x[i-2]-x[i-1] < xtol ) ) { 
      return; 
    }

    // If the spacing isn't fine enough, we'll bisect the grid
    //
    // x = [   mu1      mu2          mu3            0    0   0 ... ]
    // y = [ s(mu1)   s(mu2)       s(mu3)           0    0   0 ... ]
    //
    //           will be turned into
    //
    // x = [   mu1      mu2      0.5*(mu2+mu3)     mu3   0   0 ... ]
    // y = [ s(mu1)   s(mu2)   s(0.5*(mu2+mu3))  s(mu3)  0   0 ... ]
    //
    // where s(x) is the incoherent cross section [Eq. 225] evaluated at 
    // cosine x
    
    shiftOver( i, x, y, xm, yt );

  }  // do 110 inner. This corresponds to    // go to 110 
} 



auto do_110_120_130( int& i, std::vector<double>& x, std::vector<double>& y, 
  const double& e, const double& ep, const double& tev, const double& tevz, 
  const std::vector<double>& alpha, const std::vector<double>& beta,
  const std::vector<std::vector<double>>& sab, 
  const double& az, const double& az2, const int& lasym, const double& teff, 
  const double& teff2, const int& lat, const double& cliq, const double& sb, 
  const double& sb2, const int& iinc, const int& nl, 
  const double& sigmin, std::vector<double>& s, int& nbin, double& fract, 
  double& xl, int& j, double& ymax, const double& eps, const double& seep, 
  double& yl, const double& s1bb, const double& tol, 
  const double& xtol ){
  /* For this, we fill up mu values into x, and S(a,b,mu) values into y (with
   * a and b being fixed, and mu corresponding to the values in x). The grid
   * is chosen adaptively. So let's consider x = [2, -8, -16, 0, 0, 0, ... ]
   * (in reality 2-->1 and -16 --> -1, but this is just good for discussion now)
   */



  double gral, sum = 0;

  while (true){

    // Fills up x, y with mu and S(a,b,mu) values (respectively) so that they
    // are close enough to be reasonably interpolated.
    // x = [   mu1      mu2      mu3     ...    mu_i     0    0   0 ... ]
    // y = [ s(mu1)   s(mu2)   s(mu3)    ...  s(mu_i)    0    0   0 ... ]

        //if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "SUM         " << i << "     "<< (y[i-1]+yl) << "       " << (x[i-1]-xl) << std::endl;
    do_110(i, x, y, e, ep, tev, alpha, beta, sab, az, tevz, lasym, az2, 
        teff2, lat, cliq, sb, sb2, teff, iinc, xtol, tol, ymax);
        //if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "SUM         " << i << std::endl;

    // When do_100 returns, we x and y both have i-many nonzero entries
    // On the first iteration, xl = -1, and yl = S(a,b,mu=-1)

    while ( true ){

      // std::cout << 120 << std::endl;
      
      // Sum is meant to be the integral of cross section across all values of
      // mu. 
      // x = [   -1  ..   mu_j  ..   1   0  0  0 ... ]
      // y = [ s(-1) .. s(mu_j) .. s(1)  0  0  0 ... ]
      
      sum = sum + 0.5*(y[i-1]+yl)*(x[i-1]-xl);
      // If the tolerance between my i-th point and my (i-1)th point is 
      // reasonably low, then go left to check the other side.
      xl  = x[i-1];
      yl  = y[i-1];
      i   = i - 1;

      if ( i > 1 ){ break; }
      // if (i == 1) go to 120
      
      if (i != 1) { // don't go to 120 - either go to 130 or return 
        
        // if (sum > sigmin) go to 130
        if (sum > sigmin) { 
          s[0] = sum; 
          nbin  = nl - 1;
          fract = sum/nbin;
          sum   = 0;
          i     = 3;
          j     = 0;
          xl    = -1;
          gral  = 0;
          for ( int il = 1; il < nl; ++il ){ s[il] = 0; } 
          return std::tuple<double,double,bool> { gral, sum, false };
        } 
        else { 
          s[0] = 0; 
          for ( int il = 1; il < nl; ++il ){ s[il] = 0; } 
          return std::tuple<double,double,bool> { gral, sum, true };
        }


      } // don't go to 120
    } // go to 120
  } // go to 110

}

