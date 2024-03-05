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

void Stack::init() { // overload function when length not specified 
    // NOTE: Caller takes ownership of "s". 
    int length = 5 ; 
    count = 0 ; 
    LEN = length ; // set the LEN variable that is used in the whole Stack class to the length that is initialized. 
    s = new double[length] ; // set size of stack and create new variable 
    std::cout << "initialized by constructor" << std::endl ; // after the constructor was added, init() is now private, so there can only be self initialization 
    } 

void Stack::close() { 
    // NOTE: closes allocated memory (stack s) from init()
    std::cout << "closed by destructor" << std::endl ; 
    delete[] s ; // This line is done once if you only make one instance. If you copy a variable it always first create a new s with copy constructor. At the end it goes to this function via the destructor, where it can always delete s because it is always made when copying. 
    } 

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

Stack::Stack(const Stack& s_old) { // seems to work
    init(s_old.LEN) ; // use init constructor to set new s and LEN --> clean programming 3.3j
    count = s_old.count ; //still needs to set count to old count 
    // LEN = s_old.LEN ; // This is not needed because of init. 
    // s = new double[LEN] ; // This is not needed because of init. What happens if you do this instead of init --> After calling "Stack(const Stack& s_old)" all the variables LEN, count, s are new made (and closed in the end). However the pointer needs to be treated by the new-operator. This is because if you let the pointer point to the same address it points to the same memory. This way, if you change in the new class values of the stack this also changes with the old one, because essentially it is the same memory. 
    for (int i=0; i<LEN; i++) { 
        *(s+i) = *(s_old.s+i) ;
        cout << "copy: " << *(s_old.s+i) << " " << *(s+i) << endl ;
    } 
}
