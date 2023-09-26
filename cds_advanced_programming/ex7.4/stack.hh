#ifndef STACK
#define STACK
#include <deque>
using namespace std ; 

template<class T>
class Stack {

  public:
    Stack() {cout << "constructor" << endl ;} 
    int nitems() { return s.size() ; } 
    bool empty() { return (s.size()==0) ; }
    void push(T c) ;
    T pop() ;
    void inspect() ;

  private:  
    deque<T> s ;
} ; 

template<class T> // needs to be in same file as the template class (so in the same header file)
void Stack<T>::push(T c) { 
    s.push_back(c) ;
}

template<class T>
T Stack<T>::pop() { 
    if (empty()) {
        cout << "Stack::pop() Error: stack is empty" << endl ;
        return 0 ;
    }    
    T value = s.back() ;
    s.pop_back() ;
    return value ; 
}

template<class T>
void Stack<T>::inspect() { // prints the positions and values in Stack 
    if (s.size() == 0) {
        cout << "Stack is empty" << endl ;
    } else {
        typename deque<T>::iterator iter ;
        for (iter=s.begin(); iter!=s.end(); ++iter) {
            cout << *iter << endl ;
        }
    }
}


#endif

