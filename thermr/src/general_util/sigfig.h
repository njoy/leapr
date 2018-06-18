
#ifndef THERMR_SIGFIG_HH
#define THERMR_SIGFIG_HH
#include <cmath>

inline auto sigfig(double x, int ndig, int idig){
  /*--------------------------------------------------------------------
   * Adjust x to have ndig signficant figures.  If idig is not zero,
   * shade x up or down by idig in the last significant figure.
   *--------------------------------------------------------------------
   */
   double xx,aa;
   int ipwr,ii;
   double bias=1.0000000000001e0;
   xx=0;
   if (x != 0) {
      aa=log10(std::abs(x));
      ipwr=int(aa);
      if (aa < 0) ipwr=ipwr-1;
      ipwr=ndig-1-ipwr;
      ii=round(x*std::pow(10,ipwr)+std::pow(10.0,(ndig-11)));
      if (ii >= std::pow(10,ndig)) {
         ii=ii/10;
         ipwr=ipwr-1;
      }
      ii=ii+idig;
      xx=ii*std::pow(10,(-ipwr));
   }
   xx=xx*bias;
   return xx;
}

#endif
