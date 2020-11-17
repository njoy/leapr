#include <vector>
#include <cmath>

template <typename Range>
void sbfill(Range& sb, int nbt, const double& delta, const double& be, 
  Range& s, const Range& beta, int ndmax){

  double bmin = -be - (nbt-1) * delta,
         bmax = -be + (nbt-1) * delta + delta * 0.01;

  if ( 1 + (bmax-bmin) / delta > ndmax){ throw std::exception(); }
  
  int i = 0;
  size_t j = beta.size()-1;
  double current, toLeft, bet = bmin;
  bool foundRange = false, indexInRange = false; 
  
  while (bet < bmax){
    double b = std::fabs(bet);

    foundRange = false;
    while (not foundRange){
      if ( b > beta[j] ){   
        if ( j + 1 == beta.size() ){ 
          indexInRange = b < 1.00001 * beta[j] ? true : false;
          foundRange = true;
        }
        else { j = j + 1; } 
      } // desired value is to the right
      else {              
        if ( b < beta[j-1] and j > 1 ){
          j = j - 1;
        }
        else {
          indexInRange = true;
          foundRange   = true;
        } 
      } // desired value is to the left
    } // have you narrowed down to correct range

    // Interpolate in this range
    if ( indexInRange ){

      // Don't take the log of a negative number
      current = s[j]   <= 0 ? -225 : log( s[j]   );
      toLeft  = s[j-1] <= 0 ? -225 : log( s[j-1] );

      sb[i] = current + (b-beta[j])*(toLeft-current)/(beta[j-1]-beta[j]);

      if (bet > 0) { sb[i] = sb[i] - bet; }

      if ( sb[i] > -225 ){ sb[i] = exp(sb[i]); }

    } // if I am at correct index ( desired point between j - 1 and j ) or j
    else { sb[i] = 0; }
    ++i;
    bet += delta;
  }
}











