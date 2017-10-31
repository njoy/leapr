#include <iostream>
#include <vector>
#include <cmath>

std::vector<double> convol( const std::vector<double>& t1, std::vector<double>& t_last, const double& delta, const int& npn, const int& nl){
    std::vector<double> t_next(t_last.size(),0.0);
    int i1, i2;
    double f1, f2, be, cc, ckk; 

    for ( int i = 0; i < npn; ++i ){  // i iterates through t_next
        for ( int j = 0; j < t1.size(); ++j ){  // j iterates through t1, t_last
            i1 = i + j;
            i2 = i - j;
            f1 = 0;
            be = ( j )*delta;
            if ( t1[j] > 0 ){   // continue if the jth entry of t1 is positive
                                // If it's zero or negative, ignore since
                                // it won't contribute to to t_next
                if ( (i + j - 1) <=  int(nl) ){
                    f1 = t_last[ i1 ]*exp( -be );
                }
                f2 = 0;
                if ( i2 >= 0 and i2+1 <= int(nl) ){
                    f2 = t_last[ i2 ];
                }
                else if ( i2 < 0 and int(1-i2) <= int(nl) ){
                    be = -i2 * delta;
                    f2 = t_last[-i2]*exp( -be );
                };
                cc = t1[j] * ( f1 + f2 );
                if ( j == 0 or j == t1.size()-1 ){ cc = 0.5 * cc; }
                t_next[i] += cc;
            } // if
        } // for
        t_next[i] = t_next[i] * delta;

        if (t_next[i] < 1e-30){ t_next[i] = 0; }
        
        cc = t_next[i];
        be = (i)*delta;
        cc += t_next[i]*exp(-be);
        if ( i== 0 or i == t_next.size()-1 ){
            cc = 0.5 * cc;
        }
        ckk += cc;
    } // for
    ckk = ckk * delta;
    return t_next;
}


