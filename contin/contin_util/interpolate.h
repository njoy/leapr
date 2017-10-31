
double interpolate( const std::vector<double>& y_vec, const double& delta, const double& x_desired ){
    /* This function takes in some y vector (which corresponds to some evenly spaced
     * x vector, with spacing delta). It find the linear approximation at some
     * desired point, and returns that corresponding value.
     */
    int index_left = (x_desired / delta);   // This is the index just to the left
                                            // of my desired x point

    // Check to make sure that the desired x point is within range
    if ( index_left < (y_vec.size() - 1) ){
        double x_left  = index_left * delta; // tabulated x values to the left and 
        double x_right = x_left + delta;     // right of my desired point

        return y_vec[index_left] + ( x_desired - x_left )*( y_vec[index_left + 1] - y_vec[index_left] )/(delta);
    }        

    return 0.0;
}



