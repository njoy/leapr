
template <typename F, typename A>
auto terpk( const A& ska, F& delta, F& x ){
  /* Interpolates in a table of ska(kappa) for a required kappa
   */
  unsigned int i = x / delta;
  return i < ska.size() - 1 ? ska[i] + (x-i*delta) * (ska[i+1]-ska[i]) / delta :
                              1.0;
}

