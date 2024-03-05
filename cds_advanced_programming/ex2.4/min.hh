// declarations file for min
#ifndef MIN 
#define MIN
#include "max.hh" // the first time this is read. Only the second time (when you get directed from the max file) this is not read, because MIN is already defined and you go to #endif right away. 

namespace mylib {
    int min_func(int a,int b) ;
    double min_func(double a,double b) ;
    int min_func(int integers[], int n) ;

}

#endif 

// To avoid multiple inclusion you first look if "MIN" is defined, then 
// if not you define it. Everything up to #endif will then only be defined 
// once. This ensures that you don't get directed from "min.hh" to "max.hh"
// back and forth and get in a loop. 