
#include "repoz.h"

void finda( int i, int na, std::fstream& ntape, std::vector<double>& a, 
    std::vector<double>& buf, int nbuf ){
  /*--------------------------------------------------------------------
   * Buffered sequential i/o routine.
   * Find na elements of array a from
   * core buffer and associated binary tape.
   *--------------------------------------------------------------------
   *
   * So is this basically reading in first na values from ntape, and putting
   * them in the vector a? It seems like this is all that's really happening
   *
   */
   int nl,inow,k,j;

   nl=nbuf/na;
   inow=i%nl;
   if (inow == 0) inow=nl;

   if (i == 1) repoz(ntape);
   //if (i == 1) repoz(-ntape);
   
   std::string tab = "\t", newline = "\n";
   if (inow == 1){
     for ( size_t i1 = 0; i1 < buf.size(); ++i1 ){ 
       ntape >> buf[i1];
     }
   }

   k=na*(inow-1);
   for ( size_t j = 0; j < na; ++j ){
      k=k+1;
      a[j]=buf[k-1];
    }
}


