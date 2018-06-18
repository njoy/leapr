#include "iel/iel_util/terp1.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include "calcem/calcem_util/sigl.h"
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_ep_mu_util/310.h"
#include "calcem/calcem_util/e_ep_mu_util/313.h"
#include "calcem/calcem_util/e_ep_mu_util/360.h"
#include "calcem/calcem_util/e_ep_mu_util/380.h"
#include "calcem/calcem_util/e_ep_mu_util/410.h"
#include <cmath>


auto e_ep_mu( double& teff, double& teff2, std::vector<double>& scr, double& za, 
  double& awr, int& ncds, double& emax, double& cliq, int iinc, int lat, 
  std::vector<double>& esi, std::vector<double>& xsi, const int& lasym, 
  std::vector<double>& alpha, std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, const double t, const double& tol, 
  const double& az, const double& az2, const double& sb, const double& sb2, 
  int& nnl, const int& nl, const int& jmax, const int& nne, int iprint ){

  int itemp,iold,inew,ne,nex;
  double temp;
  int nr,np,nwtab,i,nl1,ltt,loc,l,nb,nw,lt,it,ie,jbeta=0,j=0,iskip,il,k,jnz,ll,jj,ii,nll;
  double uu,uum,test2,xlast,ulast,xs,b,enow,ep,tev,tevz,ylast,u,xl,yl;
  const int ngrid=118, nlmax=65, nemax=5000, mumax=300, imax=20;
  int numIters = 0;
  std::vector<std::vector<double>> y(nlmax,std::vector<double> (imax,0.0));
  std::vector<double> ex(imax,0.0), x(imax,0.0), yt(nlmax,0.0), yy(imax,0.0),
    yu(2*nemax,0.0), ubar(ngrid,0.0), p2(ngrid,0.0), p3(ngrid,0.0), p(4,0.0),
    uj(mumax,0.0), sj(mumax,0.0);
  double u2 = 0,u2last = 0,u3 = 0,u3last = 0;
  std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
    1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 2.53e-4, 2.97e-4, 3.5e-4, 4.2e-4, 5.06e-4, 
    6.15e-4, 7.5e-4, 8.7e-4, 1.012e-3, 1.23e-3, 1.5e-3, 1.8e-3, 2.03e-3, 
    2.277e-3, 2.6e-3, 3e-3, 3.5e-3, 4.048e-3, 4.5e-3, 5e-3, 5.6e-3, 6.325e-3, 
    7.2e-3, 8.1e-3, 9.108e-3, 0.01, 0.01063, 0.0115, 0.012397, 0.0133, 0.01417, 
    0.015, 0.016192, 0.0182, 0.0199, 0.020493, 0.0215, 0.0228, 0.0253, 0.028, 
    0.030613, 0.0338, 0.0365, 0.0395, 0.042757, 0.0465, 0.05, 0.056925, 0.0625, 
    0.069, 0.075, 0.081972, 0.09, 0.096, 0.1035, 0.111573, 0.12, 0.128, 0.1355, 
    0.145728, 0.16, 0.172, 0.184437, 0.2, 0.2277, 0.2510392, 0.2705304, 
    0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 
    0.5032575, 0.56, 0.625, 0.7, 0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 
    1.7, 1.855, 2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 
    4.9, 5.35, 5.85, 6.4, 7.0, 7.65, 8.4, 9.15, 9.85, 10.0 };
  const double sabflg = -225, eps = 1.e-4, tolmin = 5.e-7, therm = 0.0253, 
    break_val = 3000, em9 = 1e-9, up = 1.1, dn = 0.9, uumin = 1e-5, 
    yumin = 2e-7, nlpmx = 10, bk =8.617385e-5;
  //int mth, mfh;
  // save nwtab,sabmin,nl,nlt,nlp,nlp1,nl1,nnl,jmax,nne
  tevz = therm;


  // compute kernel and write in special mf6 energy-angle format
  // 300 continue
  std::cout << "300" << std::endl;
  // if (iform == 1) go to 510 
  // ^^ Take care of this in the upstairs function
  // This will go to the other branch ( the E-mu-E' branch )
  ltt=5;
  //math=matdp;
  //mfh=6;
  //mth=mtref;
  tev=t*bk;
  teff  *= bk;
  teff2 *= bk;
  scr[1-1]=za;
  scr[2-1]=awr;
  scr[3-1]=0;
  scr[4-1]=ltt; // temporary flag for this format
  scr[5-1]=1;
  scr[6-1]=0;
  nw=6;
  // call contio(0,0,nscr,scr,nb,nw)
  ncds=ncds+1;
  scr[1-1]=1;
  scr[2-1]=1;
  scr[3-1]=-1; // LIP=-1 indicates incoherent inelastic data
  scr[4-1]=1; // LAW=1 for incoherent inelastic data
  scr[5-1]=1;
  scr[6-1]=2;
  scr[7-1]=2;
  scr[8-1]=2;
  scr[9-1]=1.e-5;
  scr[10-1]=1;
  scr[11-1]=emax;
  scr[12-1]=1;
  nw=12;
  //call tab1io(0,0,nscr,scr,nb,nw)
  ncds=ncds+2;
  scr[1-1]=temp;
  scr[2-1]=0;
  scr[3-1]=3; // lang=3 is a special code for equally probable cosines
  scr[4-1]=1;
  scr[5-1]=1;
  scr[6-1]=nne;
  scr[7-1]=nne;
  scr[8-1]=2;
  nw=8;
  // call tab2io(0,0,nscr,scr,nb,nwtab)
  ncds=ncds+2;
  cliq=0;

  if ( iinc != 1 and sab[0][0] > sab[0][1]){
    cliq=(sab[0][0]-sab[0][1])*alpha[0]/(beta[1]*beta[1]);
  }

  // loop over given incident energy grid.
  // 305 continue
  std::cout << 305 << std::endl;
  ie=0;

  // Loop over incident energy (green line in my drawing)
  bool loopE = true;
  while (loopE){

    // 310 continue
    do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2, p3, 
        ep, jbeta, iskip, j, beta.size(), lasym, x );
     sigl( nnl, yt.size(), enow, ep, tev, alpha, beta, sab, yt, tol, az, tevz, iinc, 
         lat, lasym, az2, teff2, cliq, sb, sb2, teff );

     for ( int il = 0; il < nl; ++il ){ y[il][0] = yt[il]; } // enddo


    bool moreBeta = true;
    while (moreBeta) {
      // set up next panel
      // 311 continue
      std::cout << 311 << std::endl;
      x[2-1]=x[1-1];

      for (int il = 0; il < nl; ++il ){
        y[il][2-1] = y[il][1-1];
      } // enddo


      // do 313
      do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip);

      std::cout << 316 << std::endl;
      ep = sigfig(ep,8,0);
      x[1-1]=ep;

      sigl( nnl, nlmax, enow, ep, tev, alpha, beta, sab, yt, tol, az, 
        tevz, iinc, lat, lasym, az2, teff2, cliq, sb, sb2, teff );

      for (int il = 0; il < nl; ++il ){
        y[il][1-1]=yt[il];
      } // enddo

      // adaptive subdivision of panel
      i=2;
      // compare linear approximation to true function
  

      bool passedTest = false;
      while (not passedTest){

        //std::cout << 330 << std::endl;
        std::cout << std::setprecision(15) << 330 << "   " << i << "   " << j << "   " << jbeta << "    " << numIters <<  "     " << scr[4] << std::endl;

        if ( i != imax ){
          if ( iskip != 1 ){
            if (0.5*(y[0][i-2]+y[0][i-1])*(x[i-2]-x[i-1]) >= tolmin) {
              double xm = sigfig(0.5*(x[i-2]+x[i-1]),8,0);

              if (xm > x[i-1] and xm < x[i-2]) {
                sigl( nnl, nlmax, enow, xm, tev, alpha, beta, sab, yt, tol, az, 
                    tevz, iinc, lat, lasym, az2, teff2, cliq, sb, sb2, teff );
                uu = 0; uum = 0;
                bool goto330 = false;
                for ( int k = 1; k <= nl; ++k ){
                  double ym = terp1(x[i-1],y[k-1][i-1],x[i-2],y[k-1][i-2],xm,2);

                  if (k > 1) { 
                    uu  += yt[k-1];
                    uum += ym;
                  }
                  test2 = (k > 1) ? tol : tol*std::abs(yt[k-1]);

                  if (std::abs(yt[k-1]-ym) > test2){
                    do410( i, x, xm, nl, y, yt, j );
                    goto330 = true;
                    break;
                  } 
                } 
                if (goto330){ continue; }
                if (not goto330){
                  std::cout << 350 << std::endl;  // 350 continue
                  if (std::abs(uu-uum) > 2*tol*std::abs(uu)+uumin){
                    // point passes.  save top point in stack and continue.
                    do410( i, x, xm, nl, y, yt, j );
                    continue;
                  }
                }
              } 
            } 
          } // iskip
          else {
            iskip = 0;
          }
        } 
      
        do360(j, jmax, xsi, x, xlast, ylast, ulast, u2last, u3last, tolmin, y, 
          p2, p3, nll, nl, p, ubar, ie, i );

        do380( i, j, nl, nll, scr, x, y, em9, xlast, ylast, jnz, ulast, u2last, 
            u3last, p );

        if (i >= 2) continue;
        jbeta=jbeta+1;

        if (jbeta <= int(beta.size())) break; // go to 311

        for ( int il = 0; il < nl; ++il ){
          y[il][i-1]=0;
        } // enddo


        // test fails.  add point to stack and continue.
        std::cout << 430 << "         " << ie << std::endl;
        numIters += 1;
        //if (numIters > 1 )return;
        std::cout << "-----------------------------------------------" << std::endl;
        ++j;
        uu = 0;
        u2 = 0;
        u3 = 0;
        xsi[ie-1]  += 0.5 * (x[i-1]-xlast) * (y[0][i-1]+ylast);

        ubar[ie-1] += 0.5 * (x[i-1]-xlast) * (uu+ulast);
        p2[ie-1]   += 0.5 * (x[i-1]-xlast) * (u2+u2last);
        p3[ie-1]   += 0.5 * (x[i-1]-xlast) * (u3+u3last);
        xsi[ie-1]   = sigfig(xsi[ie-1],9,0);


        do380(i, j, nl, nll, scr, x, y, em9, xlast, ylast, jnz, ulast, u2last, u3last, p );
        if (y[0][0] != 0.0) jnz=j;
        if (jnz < j) j=jnz+1;
 

        if (iprint == 2) {
          ubar[ie-1]=ubar[ie-1]/xsi[ie-1];
          p2[ie-1]=p2[ie-1]/xsi[ie-1];
          p3[ie-1]=p3[ie-1]/xsi[ie-1];
          //write(nsyso,'(/,1x,"incident energy =",1pe13.6,&
          //             &  5x,"cross section =",1pe13.6,&
          //             &  5x,"mubar,p2,p3 =",3(1pe12.4))')&
          //                    enow,xsi(ie),ubar(ie),p2(ie),p3(ie)
          //write(nsyso,'(/,5x,"exit energy",11x,"pdf",7x,"cosines")')
          //write(nsyso,'(  3x,"---------------",5x,"-----------",2x,88("-"))')
          ll=6;
          for ( int jj = 1; jj <= j; ++jj ){
            //write(nsyso,'(2x,1pe15.8,5x,1pe12.5,0p,8f11.6)')&
            // (scr(ll+ii),ii=1,nlp)
            // if (nl1.gt.nlp) write(nsyso,'(34x,8f11.6)')&
            // (scr(ll+ii),ii=nlp1,nl1)
            ll=ll+nl1;
          }
        }
        scr[1-1]=0;
        scr[2-1]=enow;
        scr[3-1]=0;
        scr[4-1]=0;
        scr[5-1]=(nl+1)*j;
        scr[6-1]=nl+1;
        ncds=ncds+1+(j*(nl+1)+5)/6;

        //std::cout << "--------------       " << j << "    " << scr[4] << std::endl;
        //if (numIters > 0 )return;
        //std::cout << "call listio(0,0,nscr,scr,nb,nw)" << std::endl;
        loc=1;
        //while (nb != 0){
        //  loc=loc+nw;
        //  std::cout << "call moreio(0,0,nscr,scr(loc),nb,nw)" << std::endl;
        //}
        if (ie < nne) {
          // std::cout << "go to 310" << std::endl;
          passedTest = true;
          moreBeta = false;
        }
        else {
          std::cout << "go to 610" << std::endl;
          return; 
        }

      } // passedTest (going back to 330)
    } // moreBeta (311)
  } //  loopE (310)

}
      //std::cout << std::setprecision(15) << std::setw(25) << yt[0] << std::setw(25) << yt[1] << std::setw(25) << yt[2] << std::endl;
      /*
        std::cout<<std::setprecision(18)<<std::setw(25)<<ubar[0]<<std::setw(25)<<ubar[1]<<std::setw(25)<<ubar[2]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<ubar[3]<<std::setw(25)<<ubar[4]<<std::setw(25)<<ubar[5]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<ubar[6]<<std::setw(25)<<ubar[7]<<std::setw(25)<<ubar[8]<<std::endl;
        std::cout<<std::endl; 
        std::cout<<std::setprecision(18)<<std::setw(25)<<p3[0]<<std::setw(25)<<p3[1]<<std::setw(25)<<p3[2]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<p3[3]<<std::setw(25)<<p3[4]<<std::setw(25)<<p3[5]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<p3[6]<<std::setw(25)<<p3[7]<<std::setw(25)<<p3[8]<<std::endl;
        std::cout<<std::endl; 
        std::cout<<std::setprecision(18)<<std::setw(25)<<p2[0]<<std::setw(25)<<p2[1]<<std::setw(25)<<p2[2]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<p2[3]<<std::setw(25)<<p2[4]<<std::setw(25)<<p2[5]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<p2[6]<<std::setw(25)<<p2[7]<<std::setw(25)<<p2[8]<<std::endl;
        std::cout<<std::endl; 
        std::cout<<std::setprecision(18)<<std::setw(25)<<x[0]<<std::setw(25)<<x[1]<<std::setw(25)<<x[2]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<x[3]<<std::setw(25)<<x[4]<<std::setw(25)<<x[5]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<x[6]<<std::setw(25)<<x[7]<<std::setw(25)<<x[8]<<std::endl;
        std::cout<<std::endl; 
        std::cout<<std::setprecision(18)<<std::setw(25)<<y[0][0]<<std::setw(25)<<y[0][1]<<std::setw(25)<<y[0][2]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<y[0][3]<<std::setw(25)<<y[0][4]<<std::setw(25)<<y[0][5]<<std::endl;
        std::cout<<std::setprecision(18)<<std::setw(25)<<y[0][6]<<std::setw(25)<<y[0][7]<<std::setw(25)<<y[0][8]<<std::endl;
        std::cout<<std::endl; 
        std::cout<<std::setprecision(15)<<xsi[0]<<"      "<<xsi[1]<<"      "<<xsi[2]<<std::endl;
        std::cout<<std::setprecision(15)<<xsi[3]<<"      "<<xsi[4]<<"      "<<xsi[5]<<std::endl;
        std::cout<<std::setprecision(15)<<xsi[6]<<"      "<<xsi[7]<<"      "<<xsi[8]<<std::endl;
       
        */


