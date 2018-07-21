#include "endout/endout_util/tpidio.h"
#include <unsupported/Eigen/CXX11/Tensor>

auto endout( std::vector<double> bragg, int nedge, 
    std::vector<double> scr, int isym, int ilog, 
    const std::vector<double>& alpha, const std::vector<double>& beta, 
    const double& b7, const int& nss, const double& aws,
    const double& awr, const double& spr, const std::vector<double>& tempr,
    std::vector<double>& dwpix, std::vector<double>& dwp1, int iel,
    const double& twt, /*int nsp,*/ int nsc, const int& za, const int& mat,
    const int nout, Eigen::Tensor<double,3>& ssm ){
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
  std::vector<double> z[16];
  // equivalence(t[0],z[0])
  double small=1.e-9, tiny=-999.e0, tol=0.9e-7, up=1.01e0, therm=.0253e0;
  double bk = 8.617385e-5;

  int nsh = 0; // nsh is actually a global variable in endf.f90, but idk what
               // to do, I think it just track line numbers? 
 

  // because this would otherwise be uninitialized - check this out later
  sbs = 0; nb = 0; nw = 0;

  // write header (not going to do because effectively useless)

   // compute bound scattering cross sections
   sb = spr*std::pow((1+awr)/awr,2);
   if (aws != 0 ) sbs *= std::pow((1+aws)/aws,2);

   // for mixed moderators, merge ssm results
   if (nss != 0 and b7 <= 0) {
      srat=sbs/sb;
      nscr=-10;
      //call repoz(nscr)
      for ( int k = 0; k < int(tempr.size()); ++k ){
        for ( int j = 0; j < int(alpha.size()); ++j ){
            //read(-nscr) (scr(i),i=1,nbeta);
            for ( int i = 0; i < int(beta.size()); ++i ){
               ssm(i,j,k)=srat*ssm(i,j,k)+scr[i];
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
   // write(nsyso,'(//)') This is just putting a blank line in the output file. 
   nprnt=0;
   //nsp=0;
   nsc=0;
   int math=1;
   int mfh=0;
   int mth=0;
   text=' ';
   // read(text,'(16a4,a2)') (t(i),i=1,17); 
   // No idea what this does, maybe it's not important? Yikes
   tpidio(0,nout,nprnt,nb,nw,nsc,nsh,mth,mfh,math);
   math=mat;
   mfh=1;
   mth=451;
   scr[0]=za;
   scr[1]=awr;
   scr[2]=-1; // LRP
   scr[3]=0;  // LFI
   scr[4]=0;  // NLIB
   scr[5]=0;  // NMOD
   //all contio(0,nout,nprnt,scr[0],nb,nw);
   scr[0]=0;  // ELIS
   scr[1]=0;  // STA
   scr[2]=0;  // LIS
   scr[3]=0;  // LISO
   scr[4]=0;  
   scr[5]=6;  // NFOR
   //call contio(0,nout,nprnt,scr[0],nb,nw);
   scr[0]=1;  // AWI
   scr[1]=0;  // EMAX
   scr[2]=0;  // LREL
   scr[3]=0;  
   scr[4]=12; // NSUB
   scr[5]=6;  // NVER
   //call contio(0,nout,nprnt,scr[0],nb,nw);
   scr[0]=0;  // TEMP
   scr[1]=0;  
   scr[2]=0;  // LDRV
   scr[3]=0;  
   scr[4]=0;  // NWD
   scr[5]=    // NXC
     (iel==0) ? 2 : 3; 

   n=6;
   nc=0;
   idone=0;


   /*
   do while (idone == 0)
      text='$'
      read(nsysi,*) text
      if (text == '$') {
         idone=1
      else
         nc=nc+1
         read(text,'(16a4,a2)') (t(j),j=1,17)
         do j=1,17
            scr(n+j)=z(j)
         }
         n=n+17
         if (n > mscr) call error('endout',&
           &'scratch storage exceeded for hollerith data',' ')
      }
   }
   scr[4]=17*nc
   l=1
   call hdatio(0,nout,nprnt,scr(l),nb,nw)
   */
   // This is just all to print out the text comments. If I'm working with
   // testcase09, with a final moder at the end, then this will appear
   // in the 101 1451 5 --> 101 1451 12 lines. So this is not technically
   // yet in the file 7 parts.


   /*
   do while (nb != 0)
      l=l+nw
      call moreio(0,nout,nprnt,scr(l),nb,nw)
   }
   */
   scr[0]=0;   
   scr[1]=0;
   scr[2]=1;   // MF1 
   scr[3]=451; // MT1
   scr[4]=     // NC1
          (iel==0) ? 5+nc+1 : 5+nc+2;
   
   scr[5]=0;   // MOD1

   ii=6;
   if (iel != 0) {
      scr[ii]   = 0;
      scr[ii+1] = 0;
      scr[ii+2] = 7; // MF2
      scr[ii+3] = 2; // MT2
      if (iel < 0) ncards=3+(2*int(tempr.size())+4)/6;
      if (iel > 0) ncards=3+(2*nedge+4)/6;
      if (iel > 0 and int(tempr.size()) > 1)
        ncards=ncards+(int(tempr.size())-1)*(1+(nedge+5)/6);
      scr[ii+4]=ncards;
      scr[ii+5]=0;
      ii=ii+6;
   }
   scr[ii]   = 0;
   scr[ii+1] = 0;
   scr[ii+2] = 7;      // MF2
   scr[ii+3] = 4;      // MT2
   ncards=2+(2*int(alpha.size())+4)/6;
   if (int(tempr.size()) > 1) ncards=ncards+(int(tempr.size())-1)*(1+(int(alpha.size())+5)/6);
   ncards=5+int(beta.size())*ncards;
   scr[ii+4] = ncards; // NC2
   scr[ii+5] = 0;      // MOD2
   nw=2;
   /*
   if (iel != 0) nw=3
   call dictio(0,nout,nprnt,scr[0],nb,nw)
   call asend(nout,nprnt)
   call afend(nout,nprnt)

   !--write incoherent elastic part
   if (iel < 0) {
      mfh=7
      mth=2
      scr[0]=za
      scr[1]=awr
      scr[2]=2
      scr[3]=0
      scr[4]=0
      scr[5]=0
      call contio(0,nout,nprnt,scr[0],nb,nw)
      scr[0]=sb*npr
      scr[1]=0
      scr[2]=0
      scr[3]=0
      scr[4]=1
      ndw=int(tempr.size())
      if (ndw == 1) ndw=2
      scr[5]=ndw
      scr[6]=ndw
      scr[7]=2
      do i=1,ndw
         if (i <= int(tempr.size())) {
            scr(2*i+7)=tempr(i)
            scr(2*i+8)=sigfig(dwpix(i),7,0)
         else
            scr(2*i+7)=scr(2*i+5)
            scr(2*i+8)=scr(2*i+6)
         }
      }
      call tab1io(0,nout,nprnt,scr[0],nb,nw)
      call asend(nout,nprnt)
   }

   !--write coherent elastic part
   if (iel >= 1) {
      mfh=7
      mth=2
      scr[0]=za
      scr[1]=awr
      scr[2]=1
      scr[3]=0
      scr[4]=0
      scr[5]=0
      call contio(0,nout,nprnt,scr[0],nb,nw)

      !--thin out the 1/e part at high energies.
      sum=0
      suml=0
      w=dwpix[0]
      if (nss > 0 and b7 == 0) w=(dwpix[0]+dwp1[0])/2
      do j=1,nedge
         e=bragg(1+2*j-2)
         sum=sum+exp(-4*w*e)*bragg(1+2*j-1)
         if (sum-suml > tol*sum) {
            jmax=j
            suml=sum
         }
      }

      !--now output the records for each temperature
      do i=1,int(tempr.size())
         if (i == 1) {
            scr[0]=tempr(i)
            scr[1]=0
            scr[2]=int(tempr.size())-1
            scr[3]=0
            scr[4]=1
            scr[5]=jmax
            scr[6]=jmax
            scr[7]=1
            ii=8
            jj=0
            j1=0
            sum=0
            w=dwpix(i)
            if (nss > 0 and b7 == 0) w=(dwpix(i)+dwp1(i))/2
            j=0
            do while (j < nedge)
               j=j+1
               e=bragg(1+2*j-2)
               if (j <= jmax) jj=jj+2
               scr(ii+jj-1-j1)=sigfig(e,7,0)
               scr(ii+jj-j1)=sum+exp(-4*w*e)*bragg(1+2*j-1)
               sum=scr(ii+jj-j1)
               scr(ii+jj-j1)=sigfig(scr(ii+jj-j1),7,0)
               if (j >= nedge or jj-j1 >= npage) {
                  if (ii >= 8) {
                     call tab1io(0,nout,nprnt,scr[0],nb,nw)
                     ii=0
                     j1=npage
                     idone=1
                  else
                     call moreio(0,nout,nprnt,scr[0],nb,nw)
                     j1=j1+npage
                  }
               }
            }
         else
            scr[0]=tempr(i)
            scr[1]=0
            scr[2]=2
            scr[3]=0
            scr[4]=jmax
            scr[5]=0
            ii=6
            jj=0
            j1=0
            sum=0
            w=dwpix(i)
            if (nss > 0 and b7 == 0) w=(dwpix(i)+dwp1(i))/2
            do j=1,nedge
               if (j <= jmax) jj=jj+1
               e=sigfig(bragg(1+2*jj-2),7,0)
               scr(ii+jj-j1)=sum+exp(-4*w*e)*bragg(1+2*jj-1)
               sum=scr(ii+jj-j1)
               scr(ii+jj-j1)=sigfig(scr(ii+jj-j1),7,0)
               if (j >= nedge or jj-j1 >= npage) {
                  if (ii >= 6) {
                     call listio(0,nout,nprnt,scr[0],nb,nw)
                     ii=0
                     j1=npage
                  else
                     call moreio(0,nout,nprnt,scr[0],nb,nw)
                     j1=j1+npage
                  }
               }
            }
         }
      }
      call asend(nout,nprnt)
   }
*/

   // write inelastic part
   mfh=7;
   mth=4;
   scr[0]=za;   // ZA
   scr[1]=awr;  // AWR
   scr[2]=0;   
   scr[3]=lat;  // LAT
   scr[4]=isym; // LASYM
   scr[5]=0;   
   //call contio(0,nout,nprnt,scr[0],nb,nw)
   scr[0]=0; 
   scr[1]=0;
   scr[2]=0;    // LLN
   if (ilog != 0) scr[2]=1;
   scr[3]=0;
   scr[4]=6;    // NI
   if (nss > 0) scr[4]=6*(nss+1);
   scr[5]=nss;  // NS (number of nonprincial scatterers)
   scr[6]=npr*spr;                         // B(1) = M0f0
   scr[7]=beta[beta.size()-1];             // B(2) = E/kbT
   scr[8]=awr;                             // B(3) = A0
   scr[9]=sigfig(therm*beta(nbeta),7,0);   // B(4) = Emax
   scr[10]=0;                              // B(5)
   scr[11]=npr;                            // B(6) = number principal 
   if (nss != 0) {
      scr[12]=b7;                          // B(7) 
      scr[13]=mss*sps;                     // B(8) = M1f1
      scr[14]=aws;                         // B(9) = A1
      scr[15]=0;                           // B(10)
      scr[16]=0;                           // B(11)
      scr[17]=mss;                         // B(12) = M1, number atoms 
   }

   //call listio(0,nout,nprnt,scr[0],nb,nw)
   

   scr[0]=0;     //  |-----------------------------------------------
   scr[1]=0;     //  |  All these zeros are just part of the tab2 format
   scr[2]=0;     //  |  described at the end of 7.4.1 of ENDF manual
   scr[3]=0;     //  |-----------------------------------------------
   scr[4]=1;     // NR
   nbt=nbeta;
   if (isym == 1 or isym == 3) nbt=2*nbeta-1;
   scr[5]=nbt;   // NB = # Beta Values

   scr[6]=nbt;   // NB = # Beta Values       |  These are just 
   scr[7]=4;     // Second interp. value     |  interpolation values
   //call tab2io(0,nout,nprnt,scr[0],nb,nw)

   ii=0;
   /*
   do i=1,nbt
      do nt=1,int(tempr.size())
         sc=1
         if (lat == 1) sc=therm/(bk*tempr(nt))
         if (nt == 1) {
            scr[0]=tempr(nt)                           // TEMP
            if (mod(isym,2) == 0) scr[1]=beta(i)
            if (mod(isym,2) == 1 and i < nbeta)&
              scr[1]=-beta(nbeta-i+1)
            if (mod(isym,2) == 1 and i >= nbeta)&      // ith Beta value
              scr[1]=beta(i-nbeta+1)
            be=scr[1]*sc
            scr[2]=int(tempr.size())-1                 // LT
            scr[3]=0                                   
            scr[4]=1                                   // NR
            scr[5]=nalpha                              // NP (# alpha)
            scr[6]=nalpha                              // NP (# alpha)
            scr[7]=4
            do j=1,nalpha
               scr(7+2*j)=alpha(j)
               if (isym == 0) {
                  if (ilog == 0) {
                     scr(8+2*j)=ssm(i,j,nt)*exp(-be/2)
                     if (scr(8+2*j) >= small) {
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     else
                        scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                     }
                  else
                     scr(8+2*j)=tiny
                     if (ssm(i,j,nt) > 0) {
                        scr(8+2*j)=log(ssm(i,j,nt))-be/2
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     }
                  }
               else if (isym == 1) {
                  if (i < nbeta) {
                     if (ilog == 0) {
                        scr(8+2*j)=ssm(nbeta-i+1,j,nt)*exp(be/2)
                        if (scr(8+2*j) >= small) {
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        }
                     else
                        scr(8+2*j)=tiny
                        if (ssm(nbeta-i+1,j,nt) > 0) {
                           scr(8+2*j)=log(ssm(nbeta-i+1,j,nt))+be/2
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        }
                     }
                  else
                     if (ilog == 0) {
                        scr(8+2*j)=ssp(i-nbeta+1,j,nt)*exp(be/2)
                        if (scr(8+2*j) >= small) {
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        }
                     else
                        scr(8+2*j)=tiny
                        if (ssp(i-nbeta+1,j,nt) > 0) {
                           scr(8+2*j)=log(ssp(i-nbeta+1,j,nt))+be/2
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        }
                     }
                  }
               else if (isym == 2) {
                  if (ilog == 0) {
                     scr(8+2*j)=ssm(i,j,nt)
                    if (scr(8+2*j) >= small) {
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     else
                        scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                     }
                  else
                     scr(8+2*j)=tiny
                     if (ssm(i,j,nt) > 0) {
                        scr(8+2*j)=log(ssm(i,j,nt))
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     }
                  }
               else if (isym == 3) {
                  if (i < nbeta) {
                     if (ilog == 0) {
                        scr(8+2*j)=ssm(nbeta-i+1,j,nt)
                        if (scr(8+2*j) >= small) {
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        }
                     else
                        scr(8+2*j)=tiny
                        if (ssm(nbeta-i+1,j,nt) > 0) {
                           scr(8+2*j)=log(ssm(nbeta-i+1,j,nt))
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        }
                     }
                  else
                     if (ilog == 0) {
                        scr(8+2*j)=ssp(i-nbeta+1,j,nt)
                        if (scr(8+2*j) >= small) {
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        }
                     else
                        scr(8+2*j)=tiny
                        if (ssp(-nbeta+1,j,nt) > 0) {
                           scr(8+2*j)=log(ssp(i-nbeta+1,j,nt))
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        }
                     }
                  }
               }
               if (ilog == 0 and scr(8+2*j) < smin) scr(8+2*j)=0
            }
            call tab1io(0,nout,nprnt,scr[0],nb,nw)
            ll=1+nw
            do while (nb != 0)
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
               ll=ll+nw
            }
         else
            scr[0]=tempr(nt)
            if (mod(isym,2) == 0) scr[1]=beta(i)
            if (mod(isym,2) == 1 and i < nbeta)&
              scr[1]=-beta(nbeta-i+1)
            if (mod(isym,2) == 1 and i >= nbeta)&
              scr[1]=beta(i-nbeta+1)
            be=scr[1]*sc
            scr[2]=4
            scr[3]=0
            scr[4]=nalpha
            scr[5]=0
            do j=1,nalpha
               if (isym == 0) {
                  if (ilog == 0) {
                     scr(6+j)=ssm(i,j,nt)*exp(-be/2)
                     if (scr(6+j) >= small) {
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     else
                        scr(6+j)=sigfig(scr(6+j),6,0)
                     }
                  else
                     scr(6+j)=0
                     if (ssm(i,j,nt) > 0) {
                        scr(6+j)=log(ssm(i,j,nt))-be/2
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     }
                  }
               else if (isym == 1) {
                  if (i < nbeta) {
                     if (ilog == 0) {
                        scr(6+j)=ssm(nbeta-i+1,j,nt)*exp(be/2)
                        if (scr(6+j) >= small) {
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        }
                     else
                        scr(6+j)=tiny
                        if (ssm(nbeta-i+1,j,nt) > 0) {
                           scr(6+j)=log(ssm(nbeta-i+1,j,nt))+be/2
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        }
                     }
                  else
                     if (ilog == 0) {
                        scr(6+j)=ssp(i-nbeta+1,j,nt)*exp(be/2)
                        if (scr(6+j) >= small) {
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        }
                     else
                        scr(6+j)=tiny
                        if (ssp(i-nbeta+1,j,nt) > 0) {
                           scr(6+j)=log(ssp(i-nbeta+1,j,nt))+be/2
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        }
                     }
                  }
               else if (isym == 2) {
                  if (ilog == 0) {
                     scr(6+j)=ssm(i,j,nt)
                     if (scr(6+j) >= small) {
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     else
                        scr(6+j)=sigfig(scr(6+j),6,0)
                     }
                  else
                     scr(6+j)=tiny
                     if (ssm(i,j,nt) > 0) {
                        scr(6+j)=log(ssm(i,j,nt))
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     }
                  }
               else if (isym == 3) {
                  if (i < nbeta) {
                     if (ilog == 0) {
                        scr(6+j)=ssm(nbeta-i+1,j,nt)
                        if (scr(6+j) >= small) {
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        }
                     else
                        scr(6+j)=tiny
                        if (ssm(nbeta-i+1,j,nt) > 0) {
                           scr(6+j)=log(ssm(nbeta-i+1,j,nt))
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        }
                     }
                  else
                     if (ilog == 0) {
                        scr(6+j)=ssp(i-nbeta+1,j,nt)
                        if (scr(6+j) >= small) {
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        }
                     else
                        scr(6+j)=tiny
                        if (ssp(-nbeta+1,j,nt) > 0) {
                           scr(6+j)=log(ssp(i-nbeta+1,j,nt))
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        }
                     }
                  }
               }
               if (ilog == 0 and scr(6+j) < smin) scr(6+j)=0
            }
            call listio(0,nout,nprnt,scr[0],nb,nw)
            ll=1
            do while (nb != 0)
               ll=ll+nw
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
            }
         }
      }
   }
   if (nss != 0 and b7 <= 0) {
      scr[0]=0
      scr[1]=0
      scr[2]=0
      scr[3]=0
      scr[4]=1
      ntf=int(tempr.size())
      scr[5]=ntf
      scr[6]=ntf
      scr[7]=2
      do i=1,ntf
         if (i <= int(tempr.size())) {
            scr(2*i+7)=sigfig(tempr(i),7,0)
            scr(2*i+8)=sigfig(tempf1(i),7,0)
         else
            scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
            scr(2*i+8)=scr(2*i+6)
         }
      }
      call tab1io(0,nout,nprnt,scr[0],nb,nw)
   }
   scr[0]=0
   scr[1]=0
   scr[2]=0
   scr[3]=0
   scr[4]=1
   ntf=int(tempr.size())
   scr[5]=ntf
   scr[6]=ntf
   scr[7]=2
   do i=1,ntf
      if (i <= int(tempr.size())) {
         scr(2*i+7)=sigfig(tempr(i),7,0)
         scr(2*i+8)=sigfig(tempf(i),7,0)
      else
         scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
         scr(2*i+8)=scr(2*i+6)
      }
   }
   call tab1io(0,nout,nprnt,scr[0],nb,nw)
   call asend(nout,nprnt)
   call afend(nout,nprnt)
   call amend(nout,nprnt)
   call atend(nout,nprnt)
   return
   end subroutine endout

   */
   std::cout << nscr << nedge << isym << ilog << std::endl;
  }

