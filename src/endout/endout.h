#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"



template <typename Range, typename Float> 
auto scaleDebyeWallerCoefficients( int numSecondaryScatterers, 
  int secondaryScatterType, Range& dwpix, Range& dwp1, const Range& temps, 
  const Float& awr, const Float& aws ){
  // display endf t-effective and debye-waller integral
  for (size_t i = 0; i < temps.size(); ++i){
    if (numSecondaryScatterers == 0 or secondaryScatterType > 0){
       dwpix[i] /= (awr*temps[i]*kb);
    }
    else {
       dwpix[i] /= (aws*temps[i]*kb);
       dwp1[i]  /= (awr*temps[i]*kb);
    }
  }
}


template <typename Range>//, typename Float> 
auto processCoherentElastic( const Range& bragg ){
    std::cout << bragg[0] << std::endl;
}


template <typename Range, typename Float>
auto endout( Range& sab, const Float awr, const Float& aws, const Float& spr, 
  const Float& sps, const Range& temps, int numSecondaryScatterers, 
  int secondaryScatterType, const Range& secondaryScatterVecThing, 
  const Range& alphas, const Range& betas, Range& dwpix, Range& dwp1, int iel,
  const Float& translationalWeight ){
  // compute bound scattering cross sections
  using std::pow;
  Float sb  = spr*pow(((1.0+awr)/awr),2);
  Float sbs = (aws == 0) ? 0 : sps*pow((1.0+aws)/aws,2);

  // for mixed moderators, merge ssm results
  if (numSecondaryScatterers != 0 and secondaryScatterType <= 0){
    Float srat=sbs/sb;
    for (size_t t = 0; t < temps.size(); ++t){
      for ( size_t a = 0; a < alphas.size(); ++a ){
        for ( size_t b = 0; b < betas.size(); ++b ){      
          sab[b+a*betas.size()] = srat*sab[b+a*betas.size()] 
                                + secondaryScatterVecThing[b+a*betas.size()];
        }
      }
    }
  }

  scaleDebyeWallerCoefficients( numSecondaryScatterers, secondaryScatterType, 
                                dwpix, dwp1, temps, awr, aws );


  if (iel == 0 and translationalWeight == 0.0){
    iel = -1;
    std::cout << " this has yet to be done! " << std::endl;
    // Write out the incoherent elastic part
  }
  else if (iel > 0){
    // Write out the coherent elastic part
    
  }

  return;
  std::cout << sab[0] << std::endl;

  /*




   !--write incoherent elastic part
   if (iel.lt.0) then
      ndw=ntempr
      if (ndw.eq.1) ndw=2
      scr(6)=ndw
      scr(7)=ndw
      scr(8)=2
      do i=1,ndw
         if (i.le.ntempr) then
            scr(2*i+7)=tempr(i)
            scr(2*i+8)=sigfig(dwpix(i),7,0)
         else
            scr(2*i+7)=scr(2*i+5)
            scr(2*i+8)=scr(2*i+6)
         endif
      enddo


   !--write coherent elastic part
   if (iel.ge.1) then
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
      if (nss.gt.0.and.b7.eq.zero) w=(dwpix(1)+dwp1(1))/2
      do j=1,nedge
         e=bragg(1+2*j-2)
         sum=sum+exp(-4*w*e)*bragg(1+2*j-1)
         if (sum-suml.gt.tol*sum) then
            jmax=j
            suml=sum
         endif
      enddo

      !--now output the records for each temperature
      do i=1,ntempr
         if (i.eq.1) then
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
            if (nss.gt.0.and.b7.eq.zero) w=(dwpix(i)+dwp1(i))/2
            j=0
            do while (j.lt.nedge)
               j=j+1
               e=bragg(1+2*j-2)
               if (j.le.jmax) jj=jj+2
               scr(ii+jj-1-j1)=sigfig(e,7,0)
               scr(ii+jj-j1)=sum+exp(-4*w*e)*bragg(1+2*j-1)
               sum=scr(ii+jj-j1)
               scr(ii+jj-j1)=sigfig(scr(ii+jj-j1),7,0)
               if (j.ge.nedge.or.jj-j1.ge.npage) then
                  if (ii.ge.8) then
                     call tab1io(0,nout,nprnt,scr(1),nb,nw)
                     ii=0
                     j1=npage
                     idone=1
                  else
                     call moreio(0,nout,nprnt,scr(1),nb,nw)
                     j1=j1+npage
                  endif
               endif
            enddo
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
            if (nss.gt.0.and.b7.eq.zero) w=(dwpix(i)+dwp1(i))/2
            do j=1,nedge
               if (j.le.jmax) jj=jj+1
               e=sigfig(bragg(1+2*jj-2),7,0)
               scr(ii+jj-j1)=sum+exp(-4*w*e)*bragg(1+2*jj-1)
               sum=scr(ii+jj-j1)
               scr(ii+jj-j1)=sigfig(scr(ii+jj-j1),7,0)
               if (j.ge.nedge.or.jj-j1.ge.npage) then
                  if (ii.ge.6) then
                     call listio(0,nout,nprnt,scr(1),nb,nw)
                     ii=0
                     j1=npage
                  else
                     call moreio(0,nout,nprnt,scr(1),nb,nw)
                     j1=j1+npage
                  endif
               endif
            enddo
         endif
      enddo
      call asend(nout,nprnt)
   endif

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
   if (ilog.ne.0) scr(3)=1
   scr(4)=0
   scr(5)=6
   if (nss.gt.0) scr(5)=6*(nss+1)
   scr(6)=nss
   scr(7)=npr*spr
   scr(8)=beta(nbeta)
   scr(9)=awr
   scr(10)=sigfig(therm*beta(nbeta),7,0)
   scr(11)=0
   scr(12)=npr
   if (nss.ne.0) then
      scr(13)=b7
      scr(14)=mss*sps
      scr(15)=aws
      scr(16)=0
      scr(17)=0
      scr(18)=mss
   endif
   call listio(0,nout,nprnt,scr(1),nb,nw)
   scr(1)=0
   scr(2)=0
   scr(3)=0
   scr(4)=0
   scr(5)=1
   nbt=nbeta
   if (isym.eq.1.or.isym.eq.3) nbt=2*nbeta-1
   scr(6)=nbt
   scr(7)=nbt
   scr(8)=4
   call tab2io(0,nout,nprnt,scr(1),nb,nw)
   ii=0
   write(*,*) nbt
   do i=1,nbt
      do nt=1,ntempr
         sc=1
         if (lat.eq.1) sc=therm/(bk*tempr(nt))
         if (nt.eq.1) then
            scr(1)=tempr(nt)
            if (mod(isym,2).eq.0) scr(2)=beta(i)
            if (mod(isym,2).eq.1.and.i.lt.nbeta)&
              scr(2)=-beta(nbeta-i+1)
            if (mod(isym,2).eq.1.and.i.ge.nbeta)&
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
               if (isym.eq.0) then
                  if (ilog.eq.0) then
                     scr(8+2*j)=ssm(i,j,nt)*exp(-be/2)
                     if (scr(8+2*j).ge.small) then
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     else
                        scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                     endif
                  else
                     scr(8+2*j)=tiny
                     if (ssm(i,j,nt).gt.zero) then
                        scr(8+2*j)=log(ssm(i,j,nt))-be/2
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     endif
                  endif
               else if (isym.eq.1) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssm(nbeta-i+1,j,nt)*exp(be/2)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssm(nbeta-i+1,j,nt))+be/2
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssp(i-nbeta+1,j,nt)*exp(be/2)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssp(i-nbeta+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssp(i-nbeta+1,j,nt))+be/2
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  endif
               else if (isym.eq.2) then
                  if (ilog.eq.0) then
                     scr(8+2*j)=ssm(i,j,nt)
                    if (scr(8+2*j).ge.small) then
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     else
                        scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                     endif
                  else
                     scr(8+2*j)=tiny
                     if (ssm(i,j,nt).gt.zero) then
                        scr(8+2*j)=log(ssm(i,j,nt))
                        scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                     endif
                  endif
               else if (isym.eq.3) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssm(nbeta-i+1,j,nt)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssm(nbeta-i+1,j,nt))
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(8+2*j)=ssp(i-nbeta+1,j,nt)
                        if (scr(8+2*j).ge.small) then
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        else
                           scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                        endif
                     else
                        scr(8+2*j)=tiny
                        if (ssp(-nbeta+1,j,nt).gt.zero) then
                           scr(8+2*j)=log(ssp(i-nbeta+1,j,nt))
                           scr(8+2*j)=sigfig(scr(8+2*j),7,0)
                        endif
                     endif
                  endif
               endif
               if (ilog.eq.0.and.scr(8+2*j).lt.smin) scr(8+2*j)=0
            enddo
            call tab1io(0,nout,nprnt,scr(1),nb,nw)
            ll=1+nw
            do while (nb.ne.0)
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
               ll=ll+nw
            enddo
         else
            scr(1)=tempr(nt)
            if (mod(isym,2).eq.0) scr(2)=beta(i)
            if (mod(isym,2).eq.1.and.i.lt.nbeta)&
              scr(2)=-beta(nbeta-i+1)
            if (mod(isym,2).eq.1.and.i.ge.nbeta)&
              scr(2)=beta(i-nbeta+1)
            be=scr(2)*sc
            scr(3)=4
            scr(4)=0
            scr(5)=nalpha
            scr(6)=0
            do j=1,nalpha
               if (isym.eq.0) then
                  if (ilog.eq.0) then
                     scr(6+j)=ssm(i,j,nt)*exp(-be/2)
                     if (scr(6+j).ge.small) then
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     else
                        scr(6+j)=sigfig(scr(6+j),6,0)
                     endif
                  else
                     scr(6+j)=0
                     if (ssm(i,j,nt).gt.zero) then
                        scr(6+j)=log(ssm(i,j,nt))-be/2
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     endif
                  endif
               else if (isym.eq.1) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(6+j)=ssm(nbeta-i+1,j,nt)*exp(be/2)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssm(nbeta-i+1,j,nt))+be/2
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(6+j)=ssp(i-nbeta+1,j,nt)*exp(be/2)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssp(i-nbeta+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssp(i-nbeta+1,j,nt))+be/2
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  endif
               else if (isym.eq.2) then
                  if (ilog.eq.0) then
                     scr(6+j)=ssm(i,j,nt)
                     if (scr(6+j).ge.small) then
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     else
                        scr(6+j)=sigfig(scr(6+j),6,0)
                     endif
                  else
                     scr(6+j)=tiny
                     if (ssm(i,j,nt).gt.zero) then
                        scr(6+j)=log(ssm(i,j,nt))
                        scr(6+j)=sigfig(scr(6+j),7,0)
                     endif
                  endif
               else if (isym.eq.3) then
                  if (i.lt.nbeta) then
                     if (ilog.eq.0) then
                        scr(6+j)=ssm(nbeta-i+1,j,nt)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssm(nbeta-i+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssm(nbeta-i+1,j,nt))
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  else
                     if (ilog.eq.0) then
                        scr(6+j)=ssp(i-nbeta+1,j,nt)
                        if (scr(6+j).ge.small) then
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        else
                           scr(6+j)=sigfig(scr(6+j),6,0)
                        endif
                     else
                        scr(6+j)=tiny
                        if (ssp(-nbeta+1,j,nt).gt.zero) then
                           scr(6+j)=log(ssp(i-nbeta+1,j,nt))
                           scr(6+j)=sigfig(scr(6+j),7,0)
                        endif
                     endif
                  endif
               endif
               if (ilog.eq.0.and.scr(6+j).lt.smin) scr(6+j)=0
            enddo
            call listio(0,nout,nprnt,scr(1),nb,nw)
            ll=1
            do while (nb.ne.0)
               ll=ll+nw
               call moreio(0,nout,nprnt,scr(ll),nb,nw)
            enddo
         endif
      enddo
   enddo
   if (nss.ne.0.and.b7.le.zero) then
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
         if (i.le.ntempr) then
            scr(2*i+7)=sigfig(tempr(i),7,0)
            scr(2*i+8)=sigfig(tempf1(i),7,0)
         else
            scr(2*i+7)=sigfig(up*scr(2*i+5),7,0)
            scr(2*i+8)=scr(2*i+6)
         endif
      enddo
      call tab1io(0,nout,nprnt,scr(1),nb,nw)
   endif

   */

}
