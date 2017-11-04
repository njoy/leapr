#include <iostream>
#include <vector>
#include <cmath>

auto sbfill(std::vector<double>& sb, int nbt, double delta,
   double be, std::vector<double>& s,std::vector<double>& betan,
   int nbeta, int nbe, int ndmax){

    double bmin=-be-(nbt-1)*delta;
    double bmax=-be+(nbt-1)*delta+delta/100;
    double shade = 1.00001e0;
    double slim = -225.e0;
    if ((1+int((bmax-bmin)/delta)) > ndmax){ 
        std::cout <<  "There has been an error" << std::endl;
    }
    int j=nbeta-1;
    int i=-1;
    double bet=bmin;
    double st, stm;
    double arg;
    while (bet < bmax){
        i=i+1;
        double b=abs(bet);

        // search for correct beta range
        int idone = 0;
        while (idone == 0){
            if ( b > betan[j] ){
                if ( j+1 == nbeta and b < shade * betan[j] ){
                    idone = 1;
                }
                else {
                    if ( j == nbeta - 1 ){
                        idone = 2;
                    }
                    else {
                        j = j + 1;
                    }
                }
            }
            else {
                if ( b > betan[j-1] ){
                    idone = 1;
                }
                else {
                    if (j == 1){
                        idone = 1;
                    }
                    else {
                        j = j - 1;
                    }
                }
            }
        }
        if ( idone == 1 ){
            if ( s[j] < 0 ){
                st = slim;
            }
            else {
                st = log( s[j] );
            }
            if ( s[j-1] < 0 ){
                stm = slim;
            }
            else {
                stm = log( s[j-1] );
            }
            sb[i] = st + (b-betan[j])*(stm-st)/(betan[j-1]-betan[j]);
            if (bet > 0) sb[i] = sb[i] - bet;
            arg = sb[i];
            sb[i] = 0;
            if (arg > slim) sb[i] = exp(arg);
        }
        else {
            sb[i] = 0;
        }

        bet = bet + delta;
    }
    return 0;

}



// #include <iostream>
// #include <vector>
// #include <cmath>


// auto sbfill( std::vector<double>& sb, int nbt, double delta, double be,
//              std::vector<double>& s, const std::vector<double>& betan, int nbe,
//              int ndmax ){

//     double bmin = -be-(nbt-1)*delta;
//     double bmax = -be+(nbt-1)*delta + delta / 100;
//     if ( 1 + (bmax-bmin)/delta > ndmax ){
//         std::cout << "OH NO! There is an error in trans's sbfill" << std::endl;
//     }
//     int nbeta = betan.size();
//     int j = nbeta-1;
//     int i = 0;
//     int idone = 0;
//     double bet = bmin;
//     double shade=1.00001;
//     double slim=-225.0;
//     double arg;
//     double stm, st;
//     if (nbe > 0){return 0;}
//     while (bet < bmax ){
//         i = i + 1;
//         auto b = abs(bet);
//         idone = 0;
//         while ( idone == 0 ){
//             if ( b > betan[j] ){
//                 if ( j+1 == nbeta and b < shade * betan[j] ){
//                     idone = 1;
//                 }
//                 else {
//                     if ( j+1 == nbeta ){
//                         idone = 2;
//                     }
//                     else {
//                         j = j + 1;
//                     }

//                 }
//             }
//             else {
//                 if ( b > betan[j-1] ){
//                     idone = 1;
//                 }
//                 else {
//                     if ( j+1 == 2 ){
//                         idone = 1;
//                     }
//                     else {
//                         j = j - 1;
//                     }
//                 }
//             }
//         }

//         if ( idone == 1 ){
//             if ( s[j] < 0 ){
//                 st = slim;
//             }
//             else{ 
//                 st = log(s[j]);
//             }
//             if ( s[j-1] < 0 ){
//                 stm = slim;
//             }
//             else{
//                 stm = log(s[j-1]);
//             }
//             sb[i] = st + (b-betan[j])*(stm-st)/(betan[j-1]-betan[j]);
//             if ( bet > 0 ) {
//                 sb[i] = sb[i] - bet;
//                 arg = sb[i];
//                 sb[i] = 0;
//             }
//             if ( arg > slim ){ sb[i] = exp(arg);}
//         }
//         else{ sb[i] = 0; }
//         bet = bet + delta;
//         bet = bet + 0.1*bmax;

//     }
//     std::cout << "------------------------\n\n" <<std::endl;
//     return 0;

// }
