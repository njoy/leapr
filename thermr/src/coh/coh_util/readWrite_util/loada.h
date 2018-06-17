#include "repoz.h"
#ifndef THERMR_LOADA_HH
#define THERMR_LOADA_HH

void loada( int& i, int& na, std::fstream& ntape, int nbuf, std::vector<double>& a, 
    std::vector<double>& buf ){
 /*--------------------------------------------------------------------
  * Buffered sequential i/o routine.
  * Store na elements of array a into
  * core buffer and associated binary tape.
  *--------------------------------------------------------------------
  *
  * This takes the elements of a and puts them into the vector buf, and also
  * into the file ntape.
  *
  */
  int nl,ix,inow,k;
  double x,xnow;

  // What exactly is all fo this for? Great question
  nl = nbuf / na;
  x  = 1.0 * abs(i) / nl;
  ix =       abs(i) / nl;
  xnow = (x-ix) * nl;
  inow = xnow == 0 ? nl : xnow;
  k=na*(inow-1);

  // ---  vector a   ---> buffer buf
  for ( int j = 0; j < na; ++j ){
     k=k+1;
     buf[k-1]=a[j];
  }

  // if (i == 1) call repoz(-ntape)
  if ( i == 1 ) repoz(ntape);

  // ---  buffer buf ---> file ntape
  // I know this is ugly I'll do a better job later 
  if (inow == nl or i < 0){
    std::string tab = "\t", newline = "\n";
    for ( size_t i1 = 0; i1 < buf.size(); ++i1 ){ 
      ntape << buf[i1] << tab;
      if ( i1%10 == 9 ){ ntape << newline; }
    }
  }
}
#endif
