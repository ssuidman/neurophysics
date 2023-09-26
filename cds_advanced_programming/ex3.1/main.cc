// Questions: 
// 1) is it normal to write things in one line, such as in the "stack.h"-file? 
// 2) where to put the indentation of public and private? 
// 3) why is in "stack.h" the variable "count" not initialized by "int count" instead? 
// 4) the "void close()" function that is called by the destructor does not seem necessary, because there are no memory leaks. How to remove count and when in the code is it removed?
// 5) how does this new way work with calling constructors in the beginning and is this now how things are always called? 


#include <iostream> 
#include "stack.h"
using namespace std ;


int main() {
  
    Stack s ;
    //   s.init() ; // initialize Stack, but not necessary after the constructor was made 
  
    // Write doubles into Stack
    int i ;
    for (i=0 ; i<100 ; i++) {
        cout << "pushing value " << i*i << " in stack" << endl ;
        s.push(i*i) ;
    }

    // Count doubles in fifo
    cout << s.nitems() << " values in stack" << endl ;
    s.inspect() ; // inspect the stack when it is full

    // Read doubles back from fifo
    while (!s.empty()) {
        double val = s.pop() ;
        cout << "popping value " << val << " from stack" << endl ;
        if (s.nitems()==50) {
            s.inspect() ; // inspect the stack when it contains 50 items
        }  // look at stack after values are "popped out"
    }

    s.inspect() ; // insepct stack when it is empty 
    return 0 ;
}

// 3.1d)    "init()" should be private, because you don't want to intialize functions 
//          yourself, but expect this to be done when creating an instance of a class. 
// 3.1g)    The stack is filling until it hits 80(=LEN) then it is full and you get the
//          error. Then it tells you that there are 80 values in stack and inspects the
//          elements. Then it pops out the last value that went in up to the 50'th 
//          element, where there is inspectation again and then up to 0 it pops out and 
//          there is a last inspectation. 
// 
// To run the files:
//      g++ -c stack.cc
//      g++ -c main.cc
//      g++ -o ouput main.o stack.o 
//      ./output
// 
// To switch between git versions:
//      git checkout t_ex31
//      git checkout -b b_ex31
//      (edit, git add, git commit)
//      git checkout master
//      git merge b_ex31 