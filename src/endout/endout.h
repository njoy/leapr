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
  std::vector<double> z(17);
  // equivalence(t(1),z(1))
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
   // write(nsyso,'(//)') This is just putting a blank line in the output file
   nprnt=0;
   //nsp=0;
   nsc=0;
   int math=1;
   int mfh=0;
   int mth=0;
   text=' ';
   // read(text,'(16a4,a2)') (t(i),i=1,17); 
   // No idea what this does, maybe it's not important? Yikes
   std::cout << "HERE" << std::endl;
   tpidio(0,nout,nprnt,nb,nw,nsc,nsh,mth,mfh,math);
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
   /*
   n=6
   nc=0
   idone=0
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
   scr(5)=17*nc
   l=1
   call hdatio(0,nout,nprnt,scr(l),nb,nw)
   do while (nb != 0)
      l=l+nw
      call moreio(0,nout,nprnt,scr(l),nb,nw)
   }
   scr(1)=0
   scr(2)=0
   scr(3)=1
   scr(4)=451
   scr(5)=5+nc+1
   if (iel != 0) scr(5)=scr(5)+1
   scr(6)=0
   ii=6
   if (iel != 0) {
      scr(ii+1)=0
      scr(ii+2)=0
      scr(ii+3)=7
      scr(ii+4)=2
      if (iel < 0) ncards=3+(2*ntempr+4)/6
      if (iel > 0) ncards=3+(2*nedge+4)/6
      if (iel > 0 and ntempr > 1)&
        ncards=ncards+(ntempr-1)*(1+(nedge+5)/6)
      scr(ii+5)=ncards
      scr(ii+6)=0
      ii=ii+6
   }
   scr(ii+1)=0
   scr(ii+2)=0
   scr(ii+3)=7
   scr(ii+4)=4
   ncards=2+(2*nalpha+4)/6
   if (ntempr > 1) ncards=ncards+(ntempr-1)*(1+(nalpha+5)/6)
   ncards=5+nbeta*ncards
   scr(ii+5)=ncards
   scr(ii+6)=0
   nw=2
   if (iel != 0) nw=3
   call dictio(0,nout,nprnt,scr(1),nb,nw)
   call asend(nout,nprnt)
   call afend(nout,nprnt)

   !--write incoherent elastic part
   if (iel < 0) {
      mfh=7
      mth=2
      scr(1)=za
      scr(2)=awr
      scr(3)=2
      scr(4)=0
      scr(5)=0
      scr(6)=0
      call contio(0,nout,nprnt,scr(1),nb,nw)
      scr(1)=sb*npr
      scr(2)=0
      scr(3)=0
      scr(4)=0
      scr(5)=1
      ndw=ntempr
      if (ndw == 1) ndw=2
      scr(6)=ndw
      scr(7)=ndw
      scr(8)=2
      do i=1,ndw
         if (i <= ntempr) {
            scr(2*i+7)=tempr(i)
            scr(2*i+8)=sigfig(dwpix(i),7,0)
         else
            scr(2*i+7)=scr(2*i+5)
            scr(2*i+8)=scr(2*i+6)
         }
      }
      call tab1io(0,nout,nprnt,scr(1),nb,nw)
      call asend(nout,nprnt)
   }

   !--write coherent elastic part
   if (iel >= 1) {
      mfh=7
      mth=2
      scr(1)=za
      scr(2)=awr
      scr(3)=1
      scr(4)=0
      scr(5)=0
      scr(6)=0
      call contio(0,nout,nprnt,scr(1),nb,nw)

      !--thin out the 1/e part at high energies.
      sum=0
      suml=0
      w=dwpix(1)
      if (nss > 0 and b7 == 0) w=(dwpix(1)+dwp1(1))/2
      do j=1,nedge
         e=bragg(1+2*j-2)
         sum=sum+exp(-4*w*e)*bragg(1+2*j-1)
         if (sum-suml > tol*sum) {
            jmax=j
            suml=sum
         }
      }

      !--now output the records for each temperature
      do i=1,ntempr
         if (i == 1) {
            scr(1)=tempr(i)
            scr(2)=0
            scr(3)=ntempr-1
            scr(4)=0
            scr(5)=1
            scr(6)=jmax
            scr(7)=jmax
            scr(8)=1
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
                     call tab1io(0,nout,nprnt,scr(1),nb,nw)
                     ii=0
                     j1=npage
                     idone=1
                  else
                     call moreio(0,nout,nprnt,scr(1),nb,nw)
                     j1=j1+npage
                  }
               }
            }
         else
            scr(1)=tempr(i)
            scr(2)=0
            scr(3)=2
            scr(4)=0
            scr(5)=jmax
            scr(6)=0
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
                     call listio(0,nout,nprnt,scr(1),nb,nw)
                     ii=0
                     j1=npage
                  else
                     call moreio(0,nout,nprnt,scr(1),nb,nw)
                     j1=j1+npage
                  }
               }
            }
         }
      }
      call asend(nout,nprnt)
   }

   !--write inelastic part
   mfh=7
   mth=4
   scr(1)=za
   scr(2)=awr
   scr(3)=0
   scr(4)=lat
   scr(5)=isym
   scr(6)=0
   call contio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   if (ilog != 0) scr(3)=1
   scr(4)=0
   scr(5)=6
   if (nss > 0) scr(5)=6*(nss+1)
   scr(6)=nss
   scr(7)=npr*spr
   scr(8)=beta(nbeta)
   scr(9)=awr
   scr(10)=sigfig(therm*beta(nbeta),7,0)
   scr(11)=0
   scr(12)=npr
   if (nss != 0) {
      scr(13)=b7
      scr(14)=mss*sps
      scr(15)=aws
      scr(16)=0
      scr(17)=0
      scr(18)=mss
   }
   call listio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   nbt=nbeta
   if (isym == 1 or isym == 3) nbt=2*nbeta-1
   scr(6)=nbt
   scr(7)=nbt
   scr(8)=4
   call tab2io(0,nout,nprnt,scr(1),nb,nw)
   ii=0
   do i=1,nbt
      do nt=1,ntempr
         sc=1
         if (lat == 1) sc=therm/(bk*tempr(nt))
         if (nt == 1) {
            scr(1)=tempr(nt)
            if (mod(isym,2) == 0) scr(2)=beta(i)
            if (mod(isym,2) == 1 and i < nbeta)&
              scr(2)=-beta(nbeta-i+1)
            if (mod(isym,2) == 1 and i >= nbeta)&
              scr(2)=beta(i-nbeta+1)
            be=scr(2)*sc
            scr(3)=ntempr-1
            scr(4)=0
            scr(5)=1
            scr(6)=nalpha
            scr(7)=nalpha
            scr(8)=4
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
            call tab1io(0,nout,nprnt,scr(1),nb,nw)
            ll=1+nw
            do while (nb != 0)
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
               ll=ll+nw
            }
         else
            scr(1)=tempr(nt)
            if (mod(isym,2) == 0) scr(2)=beta(i)
            if (mod(isym,2) == 1 and i < nbeta)&
              scr(2)=-beta(nbeta-i+1)
            if (mod(isym,2) == 1 and i >= nbeta)&
              scr(2)=beta(i-nbeta+1)
            be=scr(2)*sc
            scr(3)=4
            scr(4)=0
            scr(5)=nalpha
            scr(6)=0
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
            call listio(0,nout,nprnt,scr(1),nb,nw)
            ll=1
            do while (nb != 0)
               ll=ll+nw
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
            }
         }
      }
   }
   if (nss != 0 and b7 <= 0) {
      scr(1)=0
      scr(2)=0
      scr(3)=0
      scr(4)=0
      scr(5)=1
      ntf=ntempr
      scr(6)=ntf
      scr(7)=ntf
      scr(8)=2
      do i=1,ntf
         if (i <= ntempr) {
            scr(2*i+7)=sigfig(tempr(i),7,0)
            scr(2*i+8)=sigfig(tempf1(i),7,0)
         else
            scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
            scr(2*i+8)=scr(2*i+6)
         }
      }
      call tab1io(0,nout,nprnt,scr(1),nb,nw)
   }
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   ntf=ntempr
   scr(6)=ntf
   scr(7)=ntf
   scr(8)=2
   do i=1,ntf
      if (i <= ntempr) {
         scr(2*i+7)=sigfig(tempr(i),7,0)
         scr(2*i+8)=sigfig(tempf(i),7,0)
      else
         scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
         scr(2*i+8)=scr(2*i+6)
      }
   }
   call tab1io(0,nout,nprnt,scr(1),nb,nw)
   call asend(nout,nprnt)
   call afend(nout,nprnt)
   call amend(nout,nprnt)
   call atend(nout,nprnt)
   return
   end subroutine endout

   */
   std::cout << nscr << nedge << isym << ilog << std::endl;
  }

