#include <iostream>
using namespace std ;

const int LEN = 80 ; // default stack length

struct Stack {
  // Implementation
  double s[LEN] ;
  int count ;

  // Interface
  void init() { count = 0 ; }
  int nitems() { return count ; }
  bool full() { return (count==LEN) ; }
  bool empty() { return (count==0) ; }

  void push(double c) { 
    if (full()) {
      cout << "Stack::push() Error: stack is full" << endl ;
      return ;
    }
    s[count++] = c ;
  }
  
  double pop() { 
    if (empty()) {
      cout << "Stack::pop() Error: stack is empty" << endl ;
      return 0 ;
    }    
    return s[--count] ;
  }
} ;


int main() {
  
  Stack s ;
  s.init() ; // initialize Stack
  
  // Write doubles into Stack
  int i ;
  for (i=0 ; i<10 ; i++) {
    cout << "pushing value " << i*i << " in stack" << endl ;
    s.push(i*i) ;
  }
  
  // Count doubles in fifo
  cout << s.nitems() << " value in stack" << endl ;
  
  // Read doubles back from fifo
  while (!s.empty()) {
    double val = s.pop() ;
    cout << "popping value " << val << " from stack" << endl ;
  }

  return 0 ;
}
