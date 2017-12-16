


double tausq( int m1, int m2, int m3, double c1, double c2, double twopis ){
  return (c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis;
}

double taufcc( int m1, int m2, int m3, double c1, double twothd, double twopis ){
  return c1*(m1*m1+m2*m2+m3*m3+twothd*m1*m2+twothd*m1*m3-twothd*m2*m3)*twopis;
}

double taubcc( int m1, int m2, int m3, double c1, double twopis ){
  return c1*(m1*m1+m2*m2+m3*m3+m1*m2+m2*m3+m1*m3)*twopis;
}


