// definitions file for min
#include <iostream> 
#include "min.hh" // here you are directed to the "min.hh" file

namespace mylib {

    int min_func(int a,int b) { // int min function
    std::cout << "int function" << std::endl ;
    int minimum ; 
    if (a<=b) {
        minimum = a ;
    } else {
        minimum = b ;
    }
    return minimum ;
}

    double min_func(double a,double b) { // double min function
        std::cout << "min double function" << std::endl ;
        double minimum ; 
        if (a<=b) {
            minimum = a ;
        } else {
            minimum = b ;
        }
        return minimum ;
    }

    int min_func(int integers[], int n) { // integer min function
        std::cout << "min int list function" << std::endl ;
        int minimum = integers[0] ;
        for (int i=0; i<n; i++) {
            if (integers[i] < minimum) {
                minimum = integers[i] ;
            }
        }
        return minimum ; 
    }
        
}
