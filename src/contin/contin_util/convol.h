
auto convol( const std::vector<double>& t1, const std::vector<double>& t2, 
             const double& delta ){

  std::vector<double> t3(t2.size(),0.0);

  int i1, i2, len_t1 = int(t1.size()), len_t2 = int(t2.size());
  double f1, f2;

  for ( int i = 0; i < len_t2; ++i ){    // i iterates through t3
    for ( int j = 0; j < len_t1; ++j ){  // j iterates through t1, t2
      i1 = i + j;
      i2 = i - j;
      if ( t1[j] > 0 ){
        // Convolution will only be significant if t1[j] is not zero
    
        f1 = ( i1 - 1 > len_t1 ) ? 0 : t2[i1]*exp(-j*delta);

        if      ( i2 >= 0 and  i2-1 < len_t1 ){ f2 = t2[ i2];                  }
        else if ( i2 <  0 and -i2-3 < len_t1 ){ f2 = t2[-i2] * exp( i2*delta );}
        else                                  { f2 = 0;                        }

        // If at one of the endpoints, only give half contribution
        t3[i] += ( j == 0 or j == len_t1 - 1 ) ? t1[j] * ( f1 + f2 ) * 0.5 :
                                                 t1[j] * ( f1 + f2 );
      }

    } // for

    t3[i] = ( t3[i] * delta < 1e-30 ) ? 0 : t3[i] * delta;

  } // for
  return t3;
}


