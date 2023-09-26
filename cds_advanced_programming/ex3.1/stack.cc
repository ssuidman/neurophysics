#include <iostream>
#include "stack.h"
using namespace std ;



void Stack::push(double c) { 
    if (full()) { 
        cout << "Stack::push() Error: stack is full" << endl ; 
        return ; 
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

