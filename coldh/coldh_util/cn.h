#include <vector>

auto cn( int jj, int ll, int nn ){
  /* Calculates Clebsch-Gordon coefficients for cold hydrogen or 
   * deuterium calculation
   */
  int i, kdet, kdel, ka1, ka2, ka3, ka4, kb1, kb2, kb3, kb4, iwign;
  double s, fact, zi, a1, a2, a3, a4, b1, b2, b3, b4, rat, wign;
  kdet = (jj + ll + nn) / 2;
  std::cout << kdet*2  << std::endl;
  std::cout << jj + ll + nn  << std::endl;
  std::cout << "  "  << std::endl;
  kdel = jj + ll + nn - 2 * kdet;
  if (kdel == 0) {
    ka1 = jj + ll + nn;
    ka2 = jj + ll - nn;
    ka3 = jj - ll + nn;
    ka4 = ll - jj + nn;
    kb1 = ka1 / 2;
    kb2 = ka2 / 2;
    kb3 = ka3 / 2;
    kb4 = ka4 / 2;
    s = 0;
    fact = 1;
    for ( auto i = 0; i < ka1; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0) { 
      fact = exp(s); 
    }
    a1 = sqrt(fact);

    s = 0;
    fact = 1;
    for ( auto i = 0; i < ka2; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0) { 
      fact = exp(s);
    }
    a2 = sqrt(fact);

    s = 0;
    fact = 1;
    for ( auto i = 0; i < ka3; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0) { 
      fact = exp(s);
    }
    a3 = sqrt(fact);
    s = 0;
    fact = 1;
    for ( auto i = 0; i < ka4; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0) {
      fact = exp(s);
    }
    a4 = sqrt(fact);
    s = 0;
    b1 = 1;
    for ( auto i = 0; i < kb1; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0) { 
      b1 = exp(s);
    }
    s = 0;
    b2 = 1;
    for ( auto i = 0; i < kb2; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0){
      b2 = exp(s);
    }
    s = 0;
    b3 = 1;
    for ( auto i = 0; i < kb3; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0){
      b3 = exp(s);
    }

    s = 0;
    b4 = 1;
    for ( auto i = 0; i < kb4; ++i ){
      zi = i + 1;
      s = s + log(zi);
    }
    if (s > 0.0){
      b4 = exp(s);
    }

    rat = 2 * nn + 1;
    rat = rat / (jj + ll + nn + 1);
    iwign = (jj + ll - nn)/2;
    wign = std::pow(-1,iwign);
    wign = wign * sqrt(rat) * b1/a1 * a2/b2 * a3/b3 * a4/b4;
  }
  else { wign = 0; }
  return wign;
}


