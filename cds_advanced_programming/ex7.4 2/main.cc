#include <iostream> 
#include "stack.hh"
using namespace std ;


int main() {
  
    Stack<int> s ; // can call without specified length
    
    int i ; 
    for (i=0 ; i<10 ; i++) { // Write doubles into Stack
        cout << "pushing value " << i*i << " in stack" << endl ;
        s.push(i*i) ;
    } 
    cout << s.nitems() << " values in stack" << endl ; // Count doubles in fifo

    Stack<int> sclone(s) ; // clone s with new variable sclone. "sclone(s)" calls the default copy constructor. Just as when you create "Stack s" you call "Stack()", now you do "Stack sclone(s)" and call the default Stack(const Stack&). 
    s.inspect() ; // inspect the stack when it is full
    sclone.inspect() ; // inspect the copied sclone. This looks okay 

    while (!s.empty()) { // Read doubles back from fifo 
        double val = s.pop() ; // KIJKEN OF DIT IN 1 KEER KAN!!!!!!
        cout << "popping value " << val << " from stack" << endl ; 
    } 

    s.inspect() ; // insepct stack when it is empty 
    sclone.inspect() ; // inspection of sclone still looks good 

    for (i=0 ; i<5 ; i++) {
        cout << "pushing value " << i*100 << " in stack" << endl ;
        s.push(i*100) ;
    } 

    s.inspect() ; // insepct stack when it is empty 
    sclone.inspect() ; 

    return 0 ;
}


