#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"
#include "generalTools/sigfig.h"

template <typename Range, typename Float> 
auto processCoherentElastic( const Range& bragg, const Range& dwpix, 
  const Range& dwp1, int numSecondaryScatterers, int secondaryScatterType,
  int numEdges, const Float& tol, const Range& temps ){
  Float w = dwpix[0];
  if (numSecondaryScatterers > 0 and secondaryScatterType == 0){
    w = 0.5*(dwpix[0]+dwp1[0]);
  }

  Range scr (2*numEdges+10,0.0);
  std::vector<Range> totalSCR (temps.size()+1);

  int jmax = 0;
  Float e, sum  = 0.0, suml = 0.0;
  for (int j = 1; j <= numEdges; ++j){
    e = bragg[2*j-2];
    sum += exp(-4.0*w*e)*bragg[2*j-1];
    if ( sum-suml > tol*sum ){
      jmax = j;
      suml = sum;
    }
  }

  Range energies (numEdges+5, 0.0);
  //Range energies (numEdges+5, 0.0);

  sum = 0;
  for (size_t i = 0; i < temps.size(); ++i){
    if ( i == 0 ){
      Range scr (numEdges+5,0.0);
      w = dwpix[i];
      if (numSecondaryScatterers > 0 and secondaryScatterType == 0){
        w = 0.5*(dwpix[i]+dwp1[i]);
      }
      int j = 0;
      //int jj = 0;
      int jjj = 0;
      while ( j < numEdges ){
        ++j;
        e = bragg[2*j-2];
        //if ( j <= jmax ){ jj += 2; }

        if ( j <= jmax ){ jjj+= 1; }
        energies[jjj-1] = sigfig(e,7,0);

        //scr[jj-2] = sigfig(e,7,0);

        scr[jjj-1] = sum+exp(-4.0*w*e)*bragg[2*j-1];
        sum = scr[jjj-1];
        scr[jjj-1] = sigfig(scr[jjj-1],7,0);
      }
      energies.resize(jjj);
      scr.resize(jjj);
      totalSCR[i] = scr;
    }
    else {
      sum = 0.0;
      Range scr (numEdges+5,0.0);
      int jj = 0;
      w = dwpix[i];
      if (numSecondaryScatterers > 0 and secondaryScatterType == 0){
        w = 0.5*(dwpix[i]+dwp1[i]);
      }
      for (int j = 1; j <= numEdges; ++j){ 
        if (j <= jmax){ jj += 1; }
        e = sigfig(bragg[2*jj-2],7,0);
        scr[jj-1] = sum+exp(-4.0*w*e)*bragg[2*jj-1];
        sum = scr[jj-1];
        scr[jj-1] = sigfig(scr[jj-1],7,0);
      }
      scr.resize(jj);
      totalSCR[i] = scr;
    } 
  }  
  totalSCR[temps.size()] = energies;
  return totalSCR;
}

