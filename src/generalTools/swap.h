#ifndef LEAPR_GENERAL_TOOLS_SWAP
#define LEAPR_GENERAL_TOOLS_SWAP
template <typename Float>
void swap( Float& a, Float& b ){
    Float c = a;
    a = b;
    b = c;
}
#endif

