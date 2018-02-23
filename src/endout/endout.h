
auto endout( std::vector<double> bragg, int nedge, 
    std::vector<double> scr, int isym, int ilog, 
    const std::vector<double>& alpha, const std::vector<double>& beta, 
    const double& b7, const int& nss, const double& aws,
    const double& awr, const double& spr, const std::vector<double>& tempr,
    std::vector<double>& dwpix, std::vector<double>& dwp1, int iel,
    const double& twt, int nsp, int nsc, const int& za, const int& mat ){
  /*--------------------------------------------------------------------
   * ENDF output routine
   *--------------------------------------------------------------------
   */
  int mscr = scr.size();
  int maxb = bragg.size();

  double  time,sum,suml,w,e,sb,sbs,srat,sc,be;
  int i,j,k,l,nt,jmax,ii,jj,j1,ndw,nbt,ntf,ll;
  int ncards,idone,nc,n,nb,nw,nprnt,nscr;
  std::string text;
  std::string t;
  std::vector<double> z(17);
  // equivalence(t(1),z(1))
  double small=1.e-9, tiny=-999.e0, tol=0.9e-7, up=1.01e0, therm=.0253e0;
  double bk = 8.617385e-5;
 

  // write header (not going to do because effectively useless)

   // compute bound scattering cross sections
   sb = spr*std::pow((1+awr)/awr,2);
   if (aws != 0 ) sbs *= std::pow((1+aws)/aws,2);

   // for mixed moderators, merge ssm results
   if (nss != 0 and b7 <= 0) {
      srat=sbs/sb;
      nscr=-10;
      //call repoz(nscr)
      for ( size_t k = 0; k < tempr.size(); ++k ){
        for ( size_t j = 0; j < alpha.size(); ++j ){
            //read(-nscr) (scr(i),i=1,nbeta);
            for ( size_t i = 0; i < beta.size(); ++i ){
             //  ssm(i,j,k)=srat*ssm(i,j,k)+scr(i)
            }
         }
      }
   }

   // display endf t-effective and debye-waller integral
   for ( size_t i = 0; i < tempr.size(); ++i ){
      if (nss == 0 or b7 > 0) {
         dwpix[i] /= awr*tempr[i]*bk;
      }
      else {
         dwpix[i] /= aws*tempr[i]*bk;
         dwp1[i]  /= awr*tempr[i]*bk;
      }
   }

   // Writing again

   // write endf-6 file 1.
   if (iel == 0 and twt == 0) { iel=-1; }
   // write(nsyso,'(//)') This is just putting a blank line in the output file
   nprnt=0;
   nsp=0;
   nsc=0;
   int math=1;
   int mfh=0;
   int mth=0;
   text=' ';
   // read(text,'(16a4,a2)') (t(i),i=1,17); 
   // No idea what this does, maybe it's not important? Yikes
   //call tpidio(0,nout,nprnt,z,nb,nw);
   math=mat;
   mfh=1;
   mth=451;
   scr[1]=za;
   scr[2]=awr;
   scr[3]=-1;
   scr[4]=0;
   scr[5]=0;
   scr[6]=0;
   //all contio[0,nout,nprnt,scr[1],nb,nw];
   scr[1]=0;
   scr[2]=0;
   scr[3]=0;
   scr[4]=0;
   scr[5]=0;
   scr[6]=6;
   //call contio[0,nout,nprnt,scr[1],nb,nw];
   scr[1]=1;
   scr[2]=0;
   scr[3]=0;
   scr[4]=0;
   scr[5]=12;
   scr[6]=6;
   //call contio(0,nout,nprnt,scr(1],nb,nw];
   scr[1]=0;
   scr[2]=0;
   scr[3]=0;
   scr[4]=0;
   scr[5]=0;
   scr[6]=2;
   if (iel != 0) scr[6]=3;
  }

