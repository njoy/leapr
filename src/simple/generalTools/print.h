#ifndef LEAPR_SIMPLE_PRINT
#define LEAPR_SIMPLE_PRINT

void print(std::vector<double> x){
  std::cout << std::endl;
  for (auto& y : x){ std::cout << y << " "; }
  std::cout << std::endl;
}

template<typename T>
void print(T x){
  std::cout << std::endl << x << std::endl;
}
template<typename T>
void print(T x,T y){
  std::cout << std::endl << x  << "    " << y<< std::endl;
}



#endif

