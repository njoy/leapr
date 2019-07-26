#ifndef LEAPR_DISCRE_DISCREUTIL_SINT_HH
#define LEAPR_DISCRE_DISCREUTIL_SINT_HH

#include <iostream>
#include <vector>

template <typename Float, typename Range>
inline auto sint(const Float& x, const Range& bex, const Range& rdbex, 
  const Range& sex, const Range& betan, int b, const Float& alpha, 
  const Float& wt, const Float& tbart, int nbx ){
  // Interpolates in scattering function, using SCT approximation to 
  // extrapolate outside the range in memory. For discre, the point of this 
  // function is to help evaluate the S(a,b-b_k) part of Eq. 542. Note that 
  // (b-b_k) is given as input "x". If we can use SCT we will, else we'll have
  // to interpolate on log(S)
  using std::abs; 
  // Note that tbart is likely \bar{T_s}/T, which can be found via Eq. 530
  
  // Short Collision Time approximation
  // PROBLEM -- This SCT does not match Eq. 528. Please check. 
  if ( abs(x) > betan[betan.size()-1] ){
    if (alpha <= 0.0){ return 0.0; }
    Float ex = -(wt*alpha-abs(x))*(wt*alpha-abs(x))/(4*wt*alpha*tbart);
    // If x is positive, model it as if it was negative, and then multiply
    // by a e^(-x), since S(a,b) = e^-beta * S(a,-b) in Eq. 529
    if ( x > 0.0 ){ ex = ex - x; }
    return exp(ex)/(4.0*M_PI*wt*alpha*tbart); 
  } 
  
  int xLeft = 1, xMid = b+1, xRight = nbx;

  while ( xRight-xLeft > 1.0 ){
    if      (x == bex[xMid-1]){ // If desired point is in the middle of our two
      return sex[xMid-1];       // boundaries, we have our exact point
    }
    else if (x > bex[xMid-1]){  // If desired point is to the right of midpoint,
      xLeft = xMid;             // shift everything right and continue
      xMid  = (xRight-xMid)/2 + xMid;
    } 
    else {                      // If desired point is to the left of midpoint,
      xRight = xMid;            // shift everything left and continue
      xMid   = (xMid-xLeft)/2 + xLeft;  
    }

  }
  // If we've narrowed down the region to one span one grid spacing, quit
  Float ss1 = (sex[xLeft-1]  > 0.0) ? log(sex[xLeft-1] ) : -225.0;
  Float ss3 = (sex[xRight-1] > 0.0) ? log(sex[xRight-1]) : -225.0;
  Float ex = ( (bex[xRight-1]-x)*ss1 + (x-bex[xLeft-1])*ss3 ) * rdbex[xLeft-1];
  return ex <= -225.0 ? 0.0 : exp(ex);
}

#endif
