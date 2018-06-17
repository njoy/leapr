#include "calcem/calcem_util/sig.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include "calcem/calcem_util/sigl_util/beginningLoop.h"


auto do170(int& j, double& fract, double& sum, std::vector<double>& y,
  std::vector<double>& x, double& yl, double& xn, double& xl, double& f, double& disc,
  double& ytol, double& rf, int& i, double& xil, const double& sigmin){
 //std::cout << 170 << std::endl;
  j=j+1;

  double test=(fract-sum)*(y[i-1]-yl)/((x[i-1]-xl)*yl*yl);
  if (abs(test) > ytol or yl < sigmin ){
    //std::cout << 175 << std::endl;
    f=(y[i-1]-yl)*xil;
    rf=1/f;
    disc=(yl*rf)*(yl*rf)+2*(fract-sum)*rf;
    if (disc < 0.0) {
       // write(strng,'(''disc='',1p,e12.4)') disc
       // call mess('sigl',strng,'set to abs value and continue')
       disc=abs(disc);
    } // end if
    if (f > 0.0) xn=xl-(yl*rf)+sqrt(disc);
    if (f < 0.0) xn=xl-(yl*rf)-sqrt(disc);
    if (xn > xl and xn <= x[i-1]) { 
      // go to 190
      return;
    }
    if (xn > xl and xn < (x[i-1]+ytol*(x[i-1]-xl))) {
      // go to 180
      //std::cout << 180 << std::endl;
      // 180 continue
      xn=x[i-1];
      // go to 190
      return;
    } 
    std::cout << "call error('sigl','no legal solution (quadratic path).',' ')" << std::endl;

  }

  xn=xl+(fract-sum)/yl;
  if (xn > x[i-1]) {
    //std::cout << 180 << std::endl;
    xn=x[i-1];
    // go to 190
    return;
  }

  // if (xn >= xl and xn <= x(i)) go to 190
  if (xn < xl or xn > x[i-1]) {
    std::cout << "call error('sigl','no legal solution.',' ')" << std::endl;
  }
  // go to 190
  return;
}




int do250( std::vector<double>& x, std::vector<double>& y, double& xl, double& yl,
  int& i ){
   //std::cout << "250" << std::endl;
  xl=x[i-1];
  yl=y[i-1];
  i=i-1;
  if (i > 1) {
    return 150;
  }
  if (i == 1) {
    return 160;
  }

  return 260;

}


int do160(double add, std::vector<double>& x, std::vector<double>& y, double& xl,
  double& yl, int& i, double& xil, int& j, double& fract, int& nbin, double& sum,
  double& gral, double& xn, double& shade ){

  while ( true ){
    //std::cout << "160" << std::endl;
    add=0.5*(y[i-1]+yl)*(x[i-1]-xl);

    if (x[i-1] == xl) {
      // returns either 150, 160, or 260
      int what_next = do250( x, y, xl, yl, i ); 
      if ( what_next == 160 ){ continue; }
      return what_next;
    } 

    xil=1/(x[i-1]-xl);
   
    if (i == 1 and j == nbin-1) {
      //std::cout << "165" << std::endl;
      xn=x[i-1];
      j=j+1;
      return 190;

    }

    if (sum+add >= fract*shade and j < nbin-1){ 
      return 170;

       
    }
    sum=sum+add;
    gral=gral+0.5*(yl*x[i-1]-y[i-1]*xl)*(x[i-1]+xl)
      +(y[i-1]-yl)*(x[i-1]*x[i-1]+x[i-1]*xl+xl*xl)/3.0;
     
    // go to 250
    // returns either 150, 160, or 260
    int what_next = do250( x, y, xl, yl, i ); 
    if ( what_next != 160 ){ 
      // returns either 150 or 260
      return what_next;
    }
  }
}


auto sigl( int nlin, int nlmax, double e, double ep,
  double tev, std::vector<double> alpha, std::vector<double> beta,
  std::vector<std::vector<double>> sab, std::vector<double>& s, double tolin,
  double az, double tevz, int iinc, int lat, 
  int lasym, double az2, double teff2, double cliq, double sb,
  double sb2, double teff ){

 /*-------------------------------------------------------------------
  * This is called by calcem, and uses sig.
  * * * Attempted description: * * *
  * Calcem has to turn sig(E->E',mu) --> sig(E->E') by integrating 
  * over mu. So we need to subdivide the cosine range until the 
  * actual angular function (given here) is within a nice tolerance
  *-------------------------------------------------------------------
  * * * Actual description: * * *
  * Compute the cross section and legendre components or equally-
  * probable angles for the scattering from e to ep.  Uses linear
  * reconstruction of the angular distribution computed by sig.
  *-------------------------------------------------------------------
  */
  int nl,i,j,il,nbin;
  double b,seep,sum,xl,yl,ymax;
  double fract,gral,add,xil,xn,f,rf,disc,yn,xbar;
  double tol,s1bb;
  int imax=20;
  int length_p = nlin > 0 ? nlin : 0;
  std::vector<double> x(imax,0.0), y(imax,0.0), p(length_p,0.0);
  // character(60)::strng
  double xtol = 0.00001;
  double ytol = 0.001;
  double sigmin = 1.0e-32, eps = 1.0e-3;
  double shade = 0.99999999; 

  // constant factors
  
  // """ 
  // A stact is first primed with tthe point at zero, and the first point above
  // can be derived from the positive and negative values of beta from the 
  // evaluator's beta grid using Eq. 226.
  // """ 
  // (pg. 170 of manual)
  // So b is supposed to be the first nonzero value for our stack.
  b = (lat == 1 and iinc == 2) ? (ep-e)/tevz : (ep-e)/tev; // Eq. 226

  tol  = 0.5*tolin;
  nl   = abs(nlin);
  s1bb = sqrt(1+b*b);
  sum  = 0;
  i    = 3;
  xl   = -1;
  yl = sig(e,ep,xl,tev,alpha,beta,sab,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  seep = (ep == 0) ? 0.0 : 1.0/sqrt(e*ep);

  // adaptive calculation of cross section

  // Vary the cosine mu between 3 values. Calculate the corresponding
  // incoherent cross sections, and then return ymax to be the maximum
  // cross section of these three. 
  // mu1 = -1, 
  // mu2 = mu such that alpha equals sqrt(1+beta^2)
  // mu3 = 1
  ymax = adaptiveLinearization( x, y, e, ep, tev, tevz, alpha, beta, sab, az, az2, 
      lasym, teff, teff2, lat, cliq, sb, sb2, iinc, xl, eps, seep, s1bb );
  //if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )return;


  auto out = do_110_120_130( i, x, y, e, ep, tev, tevz, alpha, beta, sab,  
      az, az2, lasym, teff, teff2, lat, cliq, sb, sb2, iinc, nl, sigmin, s, 
      nbin, fract, xl, j, ymax, eps, seep, yl, s1bb, tol, xtol );

  if (std::get<2>(out)) { return; }

  gral = std::get<0>(out);
  sum  = std::get<1>(out);

  ymax = adaptiveLinearization( x, y, e, ep, tev, tevz, alpha, beta, sab, az, az2, 
      lasym, teff, teff2, lat, cliq, sb, sb2, iinc, xl, eps, seep, s1bb );


  bool go_straight_to_150_from_190 = true;
  while ( true ){ 
    if (go_straight_to_150_from_190){
      //std::cout << "150" << std::endl;
      do_110(i, x, y, e, ep, tev, alpha, beta, sab, az, tevz, lasym, az2, 
          teff2, lat, cliq, sb, sb2, teff, iinc, xtol, tol, ymax);

    }
    go_straight_to_150_from_190 = true; 

    // check bins for this panel

    int what_next = do160(add, x, y, xl, yl, i, xil, j, fract, nbin, sum, gral, 
        xn, shade );
    if (what_next == 150){ continue; }
    if (what_next == 260){ return; }
    if (what_next == 170 ){ do170(j, fract, sum, y, x, yl, xn, xl, f, disc, 
        ytol, rf, i, xil, sigmin); }

    //190 continue
    yn = yl + (y[i-1]-yl) * (xn-xl) * xil;
    gral = gral + (xn-xl)*( yl*0.5*(xn+xl) + (y[i-1]-yl)*xil*(-xl*0.5*(xn+xl) +
      (xn*xn+xn*xl+xl*xl)/3.0) );
    xbar = gral / fract;

    // compute legendre components
    if (nlin >= 0) {
       legndr(xbar,p,nl);
       for ( int il = 1; il < nl; ++il ){
          s[il]=s[il]+p[il]/nbin;
       } // end do
    }
    // output equally probable angles
    else {
       s[j]=xbar;
    } // end if

    // continue bin loop and linearization loop
    xl = xn;
    yl = yn;
    sum = 0;
    gral = 0;

    if (j == nbin) { return; }
    if (xl < x[i-1]) {
      // go to 160
      go_straight_to_150_from_190 = false; 
      continue;
    
    }
   
    // returns either 150, 160, or 260
    what_next = do250( x, y, xl, yl, i ); 
    if ( what_next == 160 ){  
      go_straight_to_150_from_190 = false; 
      continue;
    }
    if ( what_next == 150 ){
      continue;
    }
    if (what_next == 260 ){
      return;
    }

  }

}
