#ifndef STACK
#define STACK
#include "array.hh"
using namespace std ; 

template<class T>
class Stack {

  public:
    Stack(int length) : s(length,0), count(0) {cout << "constructor" << endl ;} 
    Stack() : s(5,0), count(0) {cout << "constructor" << endl ;} 
    int nitems() { return count ; } 
    bool full() { return (count==s.size()) ; }
    bool empty() { return (count==0) ; }
    void push(T c) ;
    T pop() ;
    void inspect() ;


  private:  
    Array<T> s ;
    int count ; 
} ; 


template<class T> // needs to be in same file as the template class (so in the same header file)
void Stack<T>::push(T c) { 
    if (full()) { 
        cout << "grow stack by 10" << endl ;
        s.resize(s.size()+10,0) ;
    } 
    s[count++] = c ; 
}

template<class T>
T Stack<T>::pop() { 
    if (empty()) {
        cout << "Stack::pop() Error: stack is empty" << endl ;
        return 0 ;
    }    
    return s[--count] ;
}

template<class T>
void Stack<T>::inspect() { // prints the positions and values in Stack 
    if (count == 0) {
        cout << "Stack is empty" << endl ;
    } else {
        for (int i=0; i<count; i++) {
            cout << i << " " << s[i] << endl ;
        }
    }
}



#endif

