#include <iostream> 
#include "stack.h"
using namespace std ;


int main() {
  
    Stack s(12) ; //    You can use an initial lenght of the stack,
    // Stack s ; //     but this is not necessary. 
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
    } 

    s.inspect() ; // insepct stack when it is empty 
    return 0 ;
}

// 3.2a)    You need the lenght LEN stored as a class member, because you want to change
//          the length LEN when growing. Here you need the functions inside the class to 
//          do this and you want therefore your variable from the class itself. 
// 3.2f)    It can be nice to initialize a size, because if you already know you want to 
//          allocate a big amount of data it takes time to grow the size of the buffer 
//          all the time. 
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
