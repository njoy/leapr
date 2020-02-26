
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
   using std::pow;
   using std::abs;
   xx=0;
   if (x != 0) {
      aa=log10(abs(x));
      ipwr=int(aa);
      if (aa < 0) ipwr=ipwr-1;
      ipwr=ndig-1-ipwr;
      ii=round(x*pow(10,ipwr)+pow(10.0,(ndig-11)));
      if (ii >= pow(10,ndig)) {
         ii=ii/10;
         ipwr=ipwr-1;
      }
      //std::cout << ii  << "   " << ipwr << "   " << idig << std::endl;
      ii=ii+idig;
      xx=ii*pow(10,(-ipwr));
   }
   xx=xx*bias;
   return xx;
}

#endif
