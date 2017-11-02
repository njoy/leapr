#include <iostream> 
#include <vector>
#include <cmath>


double terps( const std::vector<double>& sd, int nsd, const double& delta, 
            const double& be ){
    double yL, yR, y_val, min = -225.0;
    int xi = be/delta;
    if ( xi < 0 or xi > nsd - 1 ){ return 0.0; }
    double xL = xi * delta;
    yL = sd[xi]   <  0.0 ? min : log( sd[xi]   );
    yR = sd[xi+1] <= 0.0 ? min : log( sd[xi+1] );
    y_val = yL + ( be - xL ) * ( yR - yL ) / ( delta );
    return y_val >= min ? exp(y_val) : 0.0;
}
