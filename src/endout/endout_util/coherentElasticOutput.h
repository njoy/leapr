#ifndef ENDOUT_PROCESS_COHERENT_ELASTIC
#define ENDOUT_PROCESS_COHERENT_ELASTIC

#include <iostream>
#include "generalTools/constants.h"
#include "range/v3/all.hpp"
#include "generalTools/sigfig.h"
#include "ENDFtk.hpp"

template <typename Range, typename Float> 
auto processCoherentElastic( const Range& bragg, const Range& primaryDWF, 
  const Range& secondaryDWF, int numSecondaryScatterers, int secondaryScatterType,
  int numEdges, const Float& tol, const Range& temps ){

  bool secondScattererWithSCT = ( numSecondaryScatterers > 0 and 
                                  secondaryScatterType  == 0 );

  Float w = secondScattererWithSCT ? (primaryDWF[0]+secondaryDWF[0])*0.5 : primaryDWF[0];

  Range scr (2*numEdges+10,0.0);
  std::vector<Range> totalSCR (temps.size());

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

  sum = 0;
  for (size_t i = 0; i < temps.size(); ++i){
    sum = 0.0;
    Range scr (numEdges+5,0.0);

    w = secondScattererWithSCT ? (primaryDWF[i]+secondaryDWF[i])*0.5 : primaryDWF[i];

    int jj = 0;

    if ( i == 0 ){
      for ( int j = 1; j <= numEdges; ++j ){
        if (j <= jmax){ ++jj; }
        e = bragg[2*j-2];
        energies[jj-1] = sigfig(e,7,0);
        sum = sum+exp(-4.0*w*e)*bragg[2*j-1];
        scr[jj-1] = sigfig(sum,7,0);
      }
      energies.resize(jj);
    }
    else {
      for (int j = 1; j <= numEdges; ++j){ 
        if (j <= jmax){ ++jj; }
        e = sigfig(bragg[2*jj-2],7,0);
        sum = sum+exp(-4.0*w*e)*bragg[2*jj-1];
        scr[jj-1] = sigfig(sum,7,0);
      }
    } 
    scr.resize(jj);
    totalSCR[i] = scr;
  }  
  return std::make_tuple(energies,totalSCR);
}

#endif
