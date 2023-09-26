#include <iostream>
using namespace std ;


// int min_func(int a,int b) { // int min function
//     cout << "int function" << endl ;
//     int minimum ; 
//     if (a<=b) {
//         minimum = a ;
//     } else {
//         minimum = b ;
//     }
//     return minimum ;
// }

double min_func(double a,double b) { // double min function
    cout << "double function" << endl ;
    double minimum ; 
    if (a<=b) {
        minimum = a ;
    } else {
        minimum = b ;
    }
    return minimum ;
}



int min_func(int integers[], int n) { // int list min function
    cout << "int list function" << endl ;
    int minimum = integers[0] ;
    for (int i=0; i<n; i++) {
        if (integers[i] < minimum) {
            minimum = integers[i] ;
        }
    }
    return minimum ; 
}

int main() {
    int a = 3 ; // create integers, doubles and int list
    int b = 2 ;
    double c = 3 ;
    double d = 4 ; 
    int e[4] = {3,5,1,3} ;
    
    int minimum_int = min_func(a,b) ; // give two integers --> uses int function
    double minimum_double = min_func(c,d) ; // give doubles --> uses double function
    int minimum_int_list = min_func(e,4) ; // give list --> uses list function
    double minimum_error = min_func(d,a) ; // giving double and integer gives error, because it is not clear for the compiler whether to use int function or double function. This problem can be fixed by shutting off the int function. Then all cases still work, bc two integers case are promoted to doubles. 

    cout << minimum_int << endl ; // print all the values for the different cases
    cout << minimum_double << endl ; 
    cout << minimum_int_list << endl ; 
    cout << minimum_error << endl ; 

    return 0 ;
}



