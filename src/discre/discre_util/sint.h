#ifndef LEAPR_DISCRE_DISCREUTIL_SINT_HH
#define LEAPR_DISCRE_DISCREUTIL_SINT_HH

#include <iostream>
#include <vector>

template <typename Float> 
auto doSCT( const Float& alpha_wgt, const Float& tbart, const Float& x){
  if (alpha_wgt <= 0.0){ return 0.0; }
  Float ex = -std::pow((alpha_wgt-abs(x)),2)/(4*alpha_wgt*tbart);
  // If x is positive, model it as if it was negative, and then multiply
  // by a e^(-x), since S(a,b) = e^-beta * S(a,-b) in Eq. 529
  if ( x > 0.0 ){ ex = ex - x; }
  return exp(ex)/(4.0*M_PI*alpha_wgt*tbart); 
}
 


template <typename Float, typename Range>
auto sint(const Float& x, const Range& bex, const Range& rdbex, 
  const Range& sex, const Range& betan, int b, const Float& alpha, 
  const Float& wt, const Float& tbart, int nbx ){
  // Interpolates in scattering function, using SCT approximation to 
  // extrapolate outside the range in memory. For discre, the point of this 
  // function is to help evaluate the S(a,b-b_k) part of Eq. 542. Note that 
  // (b-b_k) is given as input "x". If we can use SCT we will, else we'll have
  // to interpolate on log(S) (Note that tbart is proz \bar{T_s}/T in Eq. 530)
  using std::abs; 
  
  // Short Collision Time approximation
  if ( abs(x) > betan[betan.size()-1] ){ return doSCT( wt*alpha, tbart, x ); }
  
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
