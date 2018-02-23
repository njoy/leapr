
auto endout( int ntempr, std::vector<double> bragg, int nedge, 
    std::vector<double> src, int isym, int ilog ){
  /*--------------------------------------------------------------------
   * ENDF output routine
   *--------------------------------------------------------------------
   */
  int mscr = src.size();
  int maxb = bragg.size();

  double  time,sum,suml,w,e,sb,sbs,srat,sc,be;
  int i,j,k,l,nt,jmax,ii,jj,j1,ndw,nbt,ntf,ll;
  int ncards,idone,nc,n,nb,nw,nprnt,nscr;
  std::string text;
  std::string t;
  std::vector<double> z(17);
  // equivalence(t(1),z(1))
  double small=1.e-9, tiny=-999.e0, tol=0.9e-7, up=1.01e0, therm=.0253e0;

  // write header (not going to do because effectively useless)

  }

