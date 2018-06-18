#include "calcem/calcem_util/sig.h"
#include "calcem/calcem_util/sigu_util/begin_sigu.h"
#include "calcem/calcem_util/sigu_util/do150.h"


auto sigu( int nemax, const double& e, const double& u, const double& tev, 
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  const std::vector<std::vector<double>>&sab, std::vector<double>&s, 
  const double& tolin, const double& az, const double& tevz, const int& iinc, 
  const int& lat, const int& lasym, /*const double& az2, const double& teff2, */
  const double& cliq, const double& sb, const double& sb2, const double& teff ){

  /*-------------------------------------------------------------------
   * Compute the secondary energy distribution scattering for cosine u.
   * Uses linear reconstruction with the cross section from function sig.
   *-------------------------------------------------------------------
   */
   int i, j, jbeta, imax=20;
   double sum, xl, yl, xm, ym, test, yt, tol, tolmin = 1.e-6, bmax = 20;
   std::vector<double> x(imax), y(imax);
   std::cout << std::setprecision(12);

   // constant factors
   tol=tolin;
   for ( int i = 0; i < 2*nemax; ++i ){
     s[i] = 0;
   }

   double root1 = (u*sqrt(e)+sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);
   double root2 = (u*sqrt(e)-sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);

   // adaptive calculation of cross section
   sum=0;
   x[1-1]=0;
   y[1-1] = sig( e, x[1-1], xl, tev, alpha, beta, sab, az, tevz, lasym, 
      /*az2, teff2,*/ lat, cliq, sb, sb2, teff, iinc );
   jbeta = -beta.size();
   if (lasym > 0) jbeta=1;
   j = 0;
   xl = 0;
   yl = 0;

   // set up next panel
   while (true){  
     // 111 continue
     std::cout << 111 << std::endl;
     x[2-1]=x[1-1]; 
     y[2-1]=y[1-1];


     do_113_116( jbeta, lat, x, y, e, tev, tevz, root1, u, alpha, beta, 
         sab, az, lasym, /*az2,*/ teff, /*teff2,*/ cliq, sb, sb2, iinc );

     //return;
     i = 2;

     std::cout << std::setprecision(15) << x[0] << "    " << x[1] << "     " << x[2] << std::endl;
     std::cout << std::setprecision(15) << y[0] << "    " << y[1] << "     " << y[2] << std::endl;

     // compare linear approximation to true function
     // 150 continue
     bool goTo150 = true;
     while (true){
       if ( i != imax and goTo150 ){


       if ( do150( i, x, y, xm, ym, yt, test, tolmin, e, u, tev, alpha, beta, 
         sab, tevz, lasym, az, /*az2, teff2,*/ lat, cliq, sb, sb2, teff, tol, 
	 iinc) ){continue;}

       //if (do150( i, x, y, xm, ym, yt, test )){ continue; }

         /*
         std::cout << std::setprecision(20) << 150 << "     " << y[0] << std::endl;
         if (i <= 3 or 0.5*(y[i-2]+y[i-1])*(x[i-2]-x[i-1]) >= tolmin) {
           xm = 0.5*(x[i-2]+x[i-1]);
           xm = sigfig(xm,8,0);
           if (xm > x[i-1] and xm < x[i-2]){
             ym=0.5*(y[i-2]+y[i-1]);
             yt = sig( e, xm, u, tev, alpha, beta, sab, az, tevz, lasym, 
                 az2, teff2, lat, cliq, sb, sb2, teff, iinc );
             test = tol*abs(yt);
      
             if (abs(yt-ym) > test) {
               // point fails
               i=i+1;
               x[i-1]=x[i-2];
               y[i-1]=y[i-2];
               x[i-2]=xm;
               y[i-2]=yt;
               // go to 150
              // return;
               continue;
             }
           }
         }
         */
       }
       goTo150 = true;
  
       // point passes
       // 160 continue
       std::cout << 160 << "    " << i << "    " << j << "      " << sum << std::endl;
       if ( j == 220 )return;
       j=j+1;
       s[2*j+1-1]=x[i-1];
       s[2*j+2-1]=y[i-1];
       if (j > 1) sum=sum+(y[i-1]+yl)*(x[i-1]-xl);
       xl=x[i-1];
       yl=y[i-1];

       // if (j >= nemax-1) go to 170
       if (j >= nemax-1 or (jbeta > 0 and beta[jbeta-1] > bmax )) {
         std::cout << j << "    " << nemax-1 << "    " << jbeta << "     " << bmax <<std::endl;
         s[1-1]=sum;
         s[2-1]=j;
         return; 
      }

      // continue bin loop and linearization loop
      i=i-1;
      //std::cout << i << "     " << jbeta << "    " << beta.size() << std::endl;

      // if (i > 1) go to 150
      if (i > 1) {
        continue;
      }
      //std::cout << "no longer doing 150" << std::endl;
      //
      //
      if ( (unsigned) jbeta > beta.size() and i == 1 ){
        jbeta = jbeta + 1;
        goTo150 = false;
        continue;
      } 

      break;
    } 

    jbeta=jbeta+1;

    // if (jbeta <= nbeta) go to 111
    if ( (unsigned) jbeta <= beta.size()) continue;
 

    // 170 continue
    std::cout << 170 << std::endl;
    s[1-1]=sum;
    s[2-1]=j;
    return; 
  }
}


