#include <range/v3/all.hpp>
#include "generalTools/constants.h"


template <typename Range, typename RangeOfRange>//, typename Float>
auto inelasticOutput( const Range& alphas, const Range& betas, const RangeOfRange& fullSAB,
  const Range& temps, int isym, int ilog, int lat ){

  int nbt = betas.size();
  if (isym == 1 or isym == 3){ nbt = 2*betas.size()-1; }

  Range outputBetas (nbt,0.0);
  Range scr(1000,0.0);
  for (size_t b = 0; b < (unsigned) nbt; ++b){
    for (size_t t = 0; t < temps.size(); ++t){
      double sc = 1.0;
      if (lat == 1) {sc = 0.0253/(kb*temps[t]); }
      if ( t == 0 ){
        if ( isym == 0 or isym == 2){ outputBetas[1] =  betas[b]; }
        else if ( b < betas.size() ){ outputBetas[1] = -betas[betas.size()-(b+1)+1-1];  }
        else                        { outputBetas[1] =  betas[betas.size()-(b+1)+1-1];  }
        double be = outputBetas[1]*sc;

        for ( size_t a = 0; a < alphas.size(); ++a ){
          if (isym == 0){
            if (ilog == 0){
              scr[2*a] = fullSAB[t][b+a*betas.size()]*exp(-be*0.5);
            }
          }
        }
      }
    }
  }

  
}

