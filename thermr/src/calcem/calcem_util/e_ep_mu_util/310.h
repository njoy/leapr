
auto do310( int& ie, double& enow, const std::vector<double>& egrid, 
  const double& temp, const double& bk, const double& break_val, 
  const double& therm, std::vector<double>& esi, std::vector<double>& xsi, 
  std::vector<double>& ubar, std::vector<double>& p2, std::vector<double>& p3, 
  double& ep, int& jbeta, int& iskip, int& j, const int nbeta, const int& lasym, 
  std::vector<double>& x ){
     //std::cout << 310 << std::endl;
     if (temp > break_val) {
       enow = egrid[0] * 
              exp( log(                             egrid[ie] / egrid[0] ) * 
                   log( (temp*bk/therm)*egrid[egrid.size()-1] / egrid[0] ) / 
                   log(                 egrid[egrid.size()-1] / egrid[0] )
                 );
     } // endif
     else {
       enow = egrid[ie];
     }
     esi[ie]  = enow;
     xsi[ie]  = 0.0;
     ubar[ie] = 0.0;
     p2[ie]   = 0.0;
     p3[ie]   = 0.0;
     x[0]     = 0.0;
     ep       = 0.0;
     j        = 0;
     iskip    = 0;
     jbeta = lasym > 0 ? 1 : -nbeta;
     ++ie;

}


