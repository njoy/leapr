
auto copys( std::vector<std::vector<std::vector<double>>>& sab ){
  /*--------------------------------------------------------------------
   * Copy sab for principal scatterer to scratch tape on unit 10.
   *--------------------------------------------------------------------
   */
   int i,j,k,nscr;

   nscr=-10;
//   call openz(nscr,1)
   for ( size_t k = 0; k < sab.size(); ++k ){
     for ( size_t j = 0; j < sab[0].size(); ++j ){
       //  write(-nscr) (sab(i,j,k),i=1,nbeta)
     }
  }
}


