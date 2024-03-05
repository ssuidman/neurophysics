#include <iostream>
#include "stack.h"
using namespace std ;


void Stack::init(int length) {
    // NOTE: Caller takes ownership of "s".
    count = 0 ; 
    LEN = length ; // set the LEN variable that is used in the whole Stack class to the length that is initialized. 
    s = new double[length] ; // set size of stack and create new variable
    std::cout << "initialized by constructor" << std::endl ; // after the constructor was added, init() is now private, so there can only be self initialization
    } 

void Stack::init() { // overload function when lenght not specified 
    // NOTE: Caller takes ownership of "s".
    int length = 5 ;
    count = 0 ; 
    LEN = length ; // set the LEN variable that is used in the whole Stack class to the length that is initialized. 
    s = new double[length] ; // set size of stack and create new variable
    std::cout << "initialized by constructor" << std::endl ; // after the constructor was added, init() is now private, so there can only be self initialization
    } 

void Stack::close() { 
    std::cout << "closed by destructor" << std::endl ; 
    delete[] s ; // delete s
    } // how to remove count and when is removed

void Stack::push(double c) { 
    int growth_factor = 10 ;
    if (full()) { 
        grow(growth_factor) ; // if this is too big there is too much extra memory allocated. If this is too small, the funcion "grow" needs to be adapted too often. 
    } 
    s[count++] = c ; 
}

double Stack::pop() { 
    if (empty()) {
        cout << "Stack::pop() Error: stack is empty" << endl ;
        return 0 ;
    }    
    return s[--count] ;
}

void Stack::inspect() { // prints the positions and values in Stack 
    if (count == 0) {
        cout << "Stack is empty" << endl ;
    } else {
        for (int i=0; i<count; i++) {
            cout << i << " " << s[i] << endl ;
        }
    }
}

void Stack::grow(int delta) {
    LEN += delta ;
    double* s_new = new double[LEN] ; // introduces new memory allocation
    cout << "grow stack by " << delta << endl ;
    for (int i = 0 ; i<(LEN-delta); i++) {
        *(s_new+i) = *(s+i) ;
    } 
    delete [] s ; // deletes old memory allocation, so memory allocation does not pile up 
    s = s_new ; 

}

