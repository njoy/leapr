#include <iostream>
#include <vector>
#include <cmath>

template <typename Range, typename Float>
void sbfill(Range& sb, int nbt, const Float& delta, const Float& be, 
  Range& s, const Range& betan, int ndmax){

  Float bmin = -be - (nbt-1) * delta;
  Float bmax = -be + (nbt-1) * delta + delta * 0.01;

  if ( 1 + (bmax-bmin) / delta > ndmax){ throw std::exception(); }
  
  int i = 0;
  size_t j = betan.size()-1;
  Float current, toLeft, bet = bmin;
  bool foundRange = false, indexInRange = false; 
  
  while (bet < bmax){
    Float b = std::abs(bet);

    foundRange = false;
    while (not foundRange){
      if ( b > betan[j] ){   
        if ( j + 1 == betan.size() ){ 
          indexInRange = b < 1.00001 * betan[j] ? true : false;
          foundRange = true;
        }
        else { j = j + 1; } 
      } // desired value is to the right
      else {              
        if ( b < betan[j-1] and j > 1 ){
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
      current = s[j]   < 0 ? -225 : log( s[j]   );
      toLeft  = s[j-1] < 0 ? -225 : log( s[j-1] );

      sb[i] = current + (b-betan[j])*(toLeft-current)/(betan[j-1]-betan[j]);

      if (bet > 0) { sb[i] = sb[i] - bet; }

      if ( sb[i] > -225 ){ sb[i] = exp(sb[i]); }

    } // if I am at correct index ( desired point between j - 1 and j ) or j
    else { sb[i] = 0; }
    ++i;
    bet += delta;
  }
}











