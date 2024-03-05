// definitions file for max
#include <iostream> 
#include "max.hh"

namespace mylib {

    double max_func(double a,double b) { // double min function
        std::cout << "max double function" << std::endl ;
        double maximum ; 
        if (a>=b) {
            maximum = a ;
        } else {
            maximum = b ;
        }
        return maximum ;
    }

    int max_func(int integers[], int n) { // integer min function
        std::cout << "max int list function" << std::endl ;
        int maximum = integers[0] ;
        for (int i=0; i<n; i++) {
            if (integers[i] > maximum) {
                maximum = integers[i] ;
            }
        }
        return maximum ; 
    }

}
