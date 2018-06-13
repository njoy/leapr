#include <fstream>
#include <iostream>

auto tpidio( int nin, int nout, int nscr, int nb, int nw, int nsc, int nsh, 
    int mth, int mfh, int math ){
  /*--------------------------------------------------------------------
   * Utility routine for ENDF bcd and binary tapes.  Read, write,
   * and/or convert the tape identification record to/from a.  If any
   * unit is zero, it is not used.  Positive units are bcd, and
   * negative units are binary.
   *--------------------------------------------------------------------
   */

  std::vector<double> rb(17);
  std::vector<double> a(17,0.0); // a is kind of weird, don't worry about it
                                 // too much, but do figure it out at some 
                                 // point
  std::ofstream noutFile;
  noutFile.open( "tape24" );
  
  std::string hb;
  //equivalence (rb(1),hb(1))
  int inin,inout,inscr,i;

  // input.
   if (nin < 0) {
      inin=std::abs(nin);
      //read(inin) math,mfh,mth,nb,nw,(a(i),i=1,17)
    } 
   else if (nin > 0) {
      //read(nin,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth
      nw=17;
      for ( int i = 0; i < 17; ++i ){
         a[i]=rb[i];
      }
   }

   // output
   nb = 0;
   inout = std::abs(nout);
   inscr = std::abs(nscr);
   if (nout < 0) {
      // write(inout) math,mfh,mth,nb,nw,(a(i),i=1,17)
   }
   if (nscr < 0) {
      // write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,17)
   }
   if (nout <= 0 and nscr <= 0) return;
   if (nout > 0) {
      for ( int i = 0; i < 17; ++i ){
         rb[i]=a[i];
      }
      //write(nout,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth,nsh
      // This writes the 1 0 0 0 at the beginning of the tape24 (nout)
      
      noutFile << std::setw(70) << math << std::setw(2) << std::right << mfh << 
                  std::setw(3)  << mth  << std::setw(5) << std::right << 
                  nsh << std::endl;

      nsh += 1;
   }
   if (nscr > 0) {
     // Not currently used 
      for ( int i = 0; i < 17; ++i ){
         rb[i]=a[i];
      }
      //write(nscr,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth,nsc
      nsc += 1;
   }
   noutFile.close();
   return;
}


