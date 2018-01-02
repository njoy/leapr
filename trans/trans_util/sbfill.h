#include <iostream>
#include <vector>
#include <cmath>

void sbfill(std::vector<double>& sb, int nbt, double delta,
  double be, std::vector<double>& s,std::vector<double>& betan, 
  int ndmax){
  
  double bmin = -be - (nbt-1) * delta;
  double bmax = -be + (nbt-1) * delta + delta * 0.01;

  if ( 1 + (bmax-bmin) / delta > ndmax){ 
    std::cout <<  "Oh no! Error in contin's sbfill." << std::endl;
    throw std::exception();
  }
  
  double slim = -225.e0;
  int i = 0, j = betan.size()-1;
  double current, toLeft, arg, bet = bmin;
  bool foundRange = false;
  bool indexInRange = false; 
  
  while (bet < bmax){
    double b=std::abs(bet);

    // search for correct beta range
    foundRange = false;
    while (not foundRange){

      // If desired point is to the right of current point, and I still have 
      // plenty of beta values to the right that I can explore, I'll just 
      // increase my index and keep searching. 
      // If I'm at last point adn my desired point is ``basically'' 
      // where I'm at, then I'll say that j is the index of my correct point 
      // (indexInRange = true). If not then I know that desired point is 
      // between j and j + 1 ( where j + 1 not valid index ). Either way if
      // I'm running out of room, foundRange gets set to true since I've 
      // narrowed down my location to either valid or invalid.
      if ( b > betan[j] ){   
        if ( j + 1 == betan.size() ){ 
          indexInRange = b < 1.00001 * betan[j] ? true : false;
          foundRange = true;
        }
        else { j = j + 1; } 
      } // desired value is to the right

      // If desired point is not to the right of current point, and also left of 
      // (j - 1)th point, then I know I can just decrease my index. If I know
      // it's between j - 1 and j, then I have it narrowed down.
      // If I know that j == 1, and I know that desired point is not to the 
      // right of j - 1 (farthest left point, 0), then I know that my point 
      // has to be at 0. So for those latter two cases, I found the correct
      // index location and I have found the correct range.
      else {              
        if ( b < betan[j-1] and j > 1 ){
          j = j - 1;
        }
        else {
          indexInRange = true;
          foundRange = true;
        } 
      } // desired value is to the left
    } // have you narrowed down to correct range

    // Interpolate in this range
    if ( indexInRange ){

      // Don't take the log of a negative number
      current = s[j] < 0 ? slim : log( s[j] );
      toLeft  = s[j-1] < 0 ? slim : log( s[j-1] );
      
      sb[i] = current + (b-betan[j])*(toLeft-current)/(betan[j-1]-betan[j]);

      if (bet > 0) { sb[i] = sb[i] - bet; }

      if ( sb[i] > slim ){ sb[i] = exp(sb[i]); }

    } // if I am at correct index ( desired point between j - 1 and j ) or j

    // If not indexInRange, sb[i] = 0. This only happens when the desired
    // value is not in the desired range.
    else { sb[i] = 0; }

    // Increase index
    i = i + 1;
    bet = bet + delta;
  }
}











