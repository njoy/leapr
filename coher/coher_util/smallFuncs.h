


double tausq( int l1, int l2, int l3, double c1, double c2 ){
  /* This is mean to sort of evaluate Eq. 558. The it does evaluate it, but
   * multiplied by an extra factor of 2pi.
   * This is the reciprocal lattice vector length for hexagonal lattice 
   * (e.g. graphite). 
   */
  return (c1*(l1*l1+l2*l2+l1*l2)+(l3*l3*c2))*4*M_PI*M_PI;
}

double taufcc( int l1, int l2, int l3, double c1 ){
	double twopis = 4*M_PI*M_PI;
	double twothd=2/3;
  return c1*(l1*l1+l2*l2+l3*l3+twothd*l1*l2+twothd*l1*l3-twothd*l2*l3)*twopis;
}

double taubcc( int l1, int l2, int l3, double c1, double twopis ){
  return c1*(l1*l1+l2*l2+l3*l3+l1*l2+l2*l3+l1*l3)*twopis;
}


