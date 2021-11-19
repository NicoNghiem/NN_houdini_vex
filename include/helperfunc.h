#ifndef __helperfunc_h__
#define __helperfunc_h_

// implement useful helper functions

function float kronecker(float i,j){
    return i==j;
}

function float compute_angle(vector x,y){
    return 2. * atan(length(x*length(y) - length(x)*y) / length(x * length(y) + length(x) * y));
}


#endif
