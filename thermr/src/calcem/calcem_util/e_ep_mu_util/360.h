
#include "coh/coh_util/sigcoh_util/legndr.h"

auto do360(int& j, const int& jmax, std::vector<double>& xsi, 
  std::vector<double>& x, const double& xlast, const double& ylast, 
  const double& ulast, const double& u2last, const double& u3last, 
  const double& tolmin, std::vector<std::vector<double>>& y, 
  std::vector<double>& p2, std::vector<double>& p3, int& nll, const int& nl, 
  std::vector<double>& p, std::vector<double>& ubar, const int ie, const int i ){

  std::cout << 360 << std::endl;
  double uu = 0, u2 = 0, u3 = 0;
  j=j+1;
  if (j >= jmax) {
    // call error('calcem','storage exceeded.',' ')" 
    throw std::exception();
  }
  if (j > 1) {
    nll = 3;
    for ( int il = 1; il < nl; ++il ){
      legndr( y[il][i-1], p, nll );
      uu += p[1];
      u2 += p[2];
      u3 += p[3];
    } // enddo
    uu *= y[0][i-1]/(nl-1);
    u2 *= y[0][i-1]/(nl-1);
    u3 *= y[0][i-1]/(nl-1);

    xsi[ie-1]  += 0.5 * (x[i-1]-xlast) * (y[0][i-1]+ylast);
    ubar[ie-1] += 0.5 * (x[i-1]-xlast) * (uu+ulast);
    p2[ie-1]   += 0.5 * (x[i-1]-xlast) * (u2+u2last);
    p3[ie-1]   += 0.5 * (x[i-1]-xlast) * (u3+u3last);
  } // endif
  if (j == 3 and xsi[ie-1] < tolmin){ j = 2; }
}



