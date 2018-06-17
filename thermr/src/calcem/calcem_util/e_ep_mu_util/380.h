
#include "coh/coh_util/sigcoh_util/legndr.h"
#include "general_util/sigfig.h"

auto do380( int& i, int& j, const int& nl, int& nll, 
  std::vector<double>& scr, std::vector<double>& x, 
  std::vector<std::vector<double>>& y, const double& em9, double& xlast, 
  double& ylast, int& jnz, double& ulast, double& u2last, double& u3last, 
  std::vector<double>& p ){
  std::cout << 380 << std::endl;
  int jscr = 7 + (j-1)*(nl+1);
  scr[jscr-1] = x[i-1];
  scr[jscr]   = (y[0][i-1]<em9) ? sigfig(y[0][i-1],8,0) : 
                                  sigfig(y[0][i-1],9,0) ;

  for ( int il = 1; il < nl; ++ il ){
    scr[il+jscr] = sigfig(y[il][i-1],9,0);
    if (std::abs(scr[il+jscr]) > 1){
      if (std::abs(scr[il+jscr]) > 1.0005 ){ 
        std::cout << "call mess('calcem',strng,'')" << std::endl;
      }
      // scr[il+jscr] /= std::abs(scr[il+jscr]); // Set to either -1 or 1
      scr[il+jscr] = scr[il+jscr] > 0 ? 1 : -1;
    }
  } // enddo

  xlast = x[i-1];
  ylast = y[0][i-1];
  if (ylast != 0) jnz=j;
  ulast  = 0;
  u2last = 0;
  u3last = 0;
  nll = 3;

  for ( int il = 1; il < nl; ++il ){
    legndr( y[il][i-1], p, nll );
    ulast  += p[1];
    u2last += p[2];
    u3last += p[3];
  } // enddo

  ulast  = ulast  * y[0][i-1]/(nl-1);
  u2last = u2last * y[0][i-1]/(nl-1);
  u3last = u3last * y[0][i-1]/(nl-1);
  i -= 1;

}



