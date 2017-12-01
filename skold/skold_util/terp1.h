#include <iostream> 
#include <vector>

auto terp1( const double& x1, const double& y1, const double& x2, 
  const double& y2, const double& x, double& y, int i ){
  /* Interpolate one point, where  (x1,y1) and (x2,y2) are the end points of 
   * the line,  (x,y) is the interpolated point, i is the interpolation code,
   * thr6 is the kinematic threshold for i = 6 (thr6 > 0)
   */

  double a, b, t;

  // make sure x2 .ne. x1
  if (x2 == x1) {
    y = y1;
    return;
  }

  // y is constant
  if (i == 1 or y2 == y1 or x == x1) {
    y = y1;
  } 

  // y is linear in x
  else if (i == 2) {
    y = y1 + (x-x1) * (y2-y1) / (x2-x1);
  }

  // y is linear in ln(x)
  else if (i == 3) {
    y = y1 + log(x/x1) * (y2-y1) / log(x2/x1);
  }

  // ln(y) is linear in x
  else if (i == 4) {
    y = y1 * exp((x-x1) * log(y2/y1) / (x2-x1));
  }

  // ln(y) is linear in ln(x)
  else if (i == 5) {
    if (y1 == 0.0) {
      y = y1;
    } else {
      y = y1 * exp(log(x/x1) * log(y2/y1) / log(x2/x1));
    }
  }

  /* This is commmented out because I have no idea what thr6 is and 
   * I doubt I will anytime soon. Shame
  // coulomb penetrability law (charged particles only)
  else if (i == 6) {
    if (y1 == 0.0) {
      y = y1;
    else
      t = sqrt(x1-thr6);
      b = log((x2*y2)/(x1*y1));
      b = b / (1/t-1/sqrt(x2-thr6));
      a = exp(b/t) * x1 * y1;
      y = (a/x) * exp(-b/sqrt(x-thr6));
    }
  }
  */
}


