#include <range/v3/all.hpp>
#include "generalTools/constants.h"
#include "generalTools/sigfig.h"
#include <iostream>


template <typename Range>
auto getSAB(const Range& fullSAB, int t, int a, int b, int nbeta){
  return fullSAB[t][b+a*nbeta];
}




template <typename Range, typename RangeOfRange >
auto getSABreadyToWrite( const RangeOfRange& fullSAB, const Range& temps, const Range& alphas, const Range& betas, int isym, int ilog, int lat, size_t b, const RangeOfRange& fullSAB_2 = {std::vector<double> (0)} ){
    std::vector<Range> toWrite (temps.size()); 
    // This is a vector of vectors where the ith entry of this is a vector of 
    // SAB (in order of increasing temperature) for the ith temperature

    auto outputBeta = betas[b];
    for (size_t t = 0; t < temps.size(); ++t){
      double sc = 1.0;
      if (lat == 1) {sc = 0.0253/(kb*temps[t]); }
      Range scr(alphas.size(),0.0);
      if ( isym == 0 or isym == 2){ outputBeta =  betas[b]; }
      else if ( b < betas.size() ){ outputBeta = -betas[betas.size()-b-1]; }
      else                        { outputBeta =  betas[b-betas.size()+1]; }
      double be = outputBeta*sc;

      for ( size_t a = 0; a < alphas.size(); ++a ){
        if (isym == 0){
          // No ncold option 
          // Symmetric scattering law
          if (ilog == 0){
            scr[a] = getSAB(fullSAB, t, a, b, betas.size())*exp(-be*0.5);
            scr[a] = (scr[a] < 1e-9) ? sigfig(scr[a],6,0)
                                     : sigfig(scr[a],7,0);
          }
          else {
            scr[a] = -999.0;
            auto sab = getSAB(fullSAB, t, a, b, betas.size());
            if ( sab > 0 ){
              scr[a] = log(sab)-be*0.5;
              scr[a] = sigfig(scr[a],7,0);
            }
          }
        } 

        if (isym == 1){
          if ( b < betas.size()-1 ){
            if (ilog == 0){
              scr[a] = getSAB(fullSAB, t, a, betas.size()-b-1, betas.size())*exp(be*0.5);
              scr[a] = (scr[a] < 1e-9) ? sigfig(scr[a],6,0)
                                       : sigfig(scr[a],7,0);
            }
            else {
              scr[a] = -999.0;
              auto sab = getSAB(fullSAB, t, a, betas.size()-b-1, betas.size());
              //std::cout << sab << std::endl;
              if ( sab > 0 ){
                scr[a] = log(sab)+be*0.5;
                scr[a] = sigfig(scr[a],7,0);
              }
            }
          }
          else {
            if (ilog == 0){
              scr[a] = getSAB(fullSAB_2, t, a, b-betas.size()+1, betas.size())*exp(be*0.5);
              scr[a] = (scr[a] < 1e-9) ? sigfig(scr[a],6,0)
                                       : sigfig(scr[a],7,0);
            }
            else {
              scr[a] = -999.0;
              auto sab = getSAB(fullSAB_2, t, a, b-betas.size()+1, betas.size());
              if ( sab > 0 ){
                scr[a] = log(sab)+be*0.5;
                scr[a] = sigfig(scr[a],7,0);
              }
            }
          }
        }



        if (isym == 2){
          if (ilog == 0){
            scr[a] = getSAB(fullSAB, t, a, b, betas.size());
            scr[a] = (scr[a] < 1e-9) ? sigfig(scr[a],6,0)
                                     : sigfig(scr[a],7,0);
          }
          else {
            scr[a] = -999.0;
            auto sab = getSAB(fullSAB, t, a, b, betas.size());
            if ( sab > 0 ){
              scr[a] = log(sab);//+be*0.5;
              scr[a] = sigfig(scr[a],7,0);
            }
          }
        }


        if (isym == 3){
          if ( b < betas.size()-1){
            if (ilog == 0){
              scr[a] = getSAB(fullSAB, t, a, betas.size()-b-1, betas.size());
              std::cout << scr[a] << std::endl;
              scr[a] = (scr[a] < 1e-9) ? sigfig(scr[a],6,0)
                                       : sigfig(scr[a],7,0);
            } 
            else {
              scr[a] = -999.0;
              auto sab = getSAB(fullSAB, t, a, betas.size()-b-1, betas.size());
              std::cout << scr[a] << std::endl;
              if ( sab > 0 ){
                scr[a] = log(sab);//+be*0.5;
                scr[a] = sigfig(scr[a],7,0);
              }
            }
          } 
          else {
            if (ilog == 0){
              scr[a] = getSAB(fullSAB_2, t, a, b-betas.size(), betas.size())*exp(be*0.5);
              scr[a] = (scr[a] < 1e-9) ? sigfig(scr[a],6,0)
                                       : sigfig(scr[a],7,0);
            } 
            else {
              scr[a] = -999.0;
              auto sab = getSAB(fullSAB_2, t, a, b-betas.size(), betas.size());
              if ( sab > 0 ){
                scr[a] = log(sab);//+be*0.5;
                scr[a] = sigfig(scr[a],7,0);
              }
            }
          }
        }


        if (ilog == 0 and scr[a] < -999.0){
          scr[a] = 0.0;
        }
      } 
      toWrite[t] = scr;
    }
  return std::make_tuple(outputBeta,toWrite); 
}







template <typename Range, typename RangeOfRange>//, typename Float>
auto inelasticOutput( const Range& alphas, const Range& betas, const RangeOfRange& fullSAB,
  const Range& temps, int isym, int ilog, int lat ){

  int nbt = betas.size();
  if (isym == 1 or isym == 3){ nbt = 2*betas.size()-1; }

  Range outputBetas (nbt,0.0);
  for (size_t b = 0; b < (unsigned) nbt; ++b){
    getSABreadyToWrite( fullSAB, temps, alphas, betas, isym, ilog, lat, b );

  }








          /*
          else if (isym == 1){
            // Ncold option 
            // Symmetric scattering law

            if (i < betas.size()){
              if (ilog == 0){
                scr[2*a] = getSAB(fullSAB, t, a, betas.size()-b, betas.size())*exp(be*0.5);
                if (scr(8+2*j) >= 1e-9) {
                  scr[2*a] = sigfig(scr[2*a],7,0);
                }
                else {
                  scr(8+2*j)=sigfig(scr(8+2*j),6,0)
                }
              }
              else {
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

          } 
          else if (isym == 2){
            // No ncold option 
            // Non-symmetric scattering law

          } 
          else if (isym == 3){
            // Ncold option 
            // Non-symmetric scattering law

          } 





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

                  */

  
}

