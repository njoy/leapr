#ifndef LEAPR_DISCRE_DISCREUTIL_SINT_HH
#define LEAPR_DISCRE_DISCREUTIL_SINT_HH

#include <iostream>
#include <vector>

inline auto sint(const double& x, const std::vector<double>& bex, 
  const std::vector<double>& rdbex, const std::vector<double>& sex,
  const std::vector<double>& betan, int b, const double& alpha,
  const double& wt, const double& tbart, int nbx ){
  // Interpolates in scattering function, using SCT approximation to 
  // extrapolate outside the range in memory. For discre, the point of this 
  // function is to help evaluate the S(a,b-b_k) part of Eq. 542. Note that 
  // (b-b_k) is given as input "x". If we can use SCT we will, else we'll have
  // to interpolate on log(S)
  
  // Note that tbart is likely \bar{T_s}/T, which can be found via Eq. 530
  
  double sv, ex, sint, pi = 3.1415926;
  
  // Short Collision Time approximation
  // PROBLEM -- This SCT does not match Eq. 528. Please check. 
 
  if ( std::abs(x) > betan[betan.size()-1] ){
    if ( alpha <= 0.0 ){ // The formula for short collision time approximation
      return 0.0;        // is only valid for positive alpha, else get complex
                         // S(a,-b) (square root)
    } 
    else {  // This implements the short collision time (SCT) approximaation
            // which is detailed in Eq. 528 of NJOY manual
      ex = -(wt*alpha-std::abs(x))*(wt*alpha-std::abs(x))/(4*wt*alpha*tbart);

      // If x is positive, model it as if it was negative, and then multiply
      // by a e^(-x), since S(a,b) = e^-beta * S(a,-b) in Eq. 529
      if ( x > 0.0 ){ ex = ex - x; }

      return exp(ex)/(4*pi*wt*alpha*tbart); // Continue with Eq. 528
    }
  } 
  
  // interpolation
  int xLeft = 1, xMid = b+1, xRight = nbx;

  // bisect for x
  int idone = 0;

  while ( idone == 0 ){
    if ( x == bex[xMid-1] ){ // If desired point is in the middle of our two
      return sex[xMid-1];    // boundaries, we have our exact point
    }
    if ( x > bex[xMid-1] ){  // If desired point is to the right of midpoint,
      xLeft = xMid;          // shift everything right and continue
      xMid = (xRight-xMid)/2 + xMid;
    } 
    else {                  // If desired point is to the left of midpoint,
      xRight = xMid;        // shift everything left and continue
      xMid = (xMid-xLeft)/2 + xLeft;  
    }

    if ( xRight-xLeft <= 1 ){ // If we've narrowed down the region to one span
      idone = 1;              // one grid spacing, quit
    }
  }
  double ss1, ss3;
  ss1 = sex[xLeft-1]  > 0.0 ? log( sex[xLeft-1]  ) : -225.0;
  ss3 = sex[xRight-1] > 0.0 ? log( sex[xRight-1] ) : -225.0;
  
  ex = ( (bex[xRight-1]-x)*ss1+(x-bex[xLeft-1])*ss3 ) * rdbex[xLeft-1];

  return ex <= -225.0 ? 0.0 : exp(ex);

}

#endif
