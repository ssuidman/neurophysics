#include <iostream> 
#include "stack.h"
using namespace std ;


int main() {
  
    Stack s ; // You can also call s(7) if you want the initial stack to be of length 7 
    // Write doubles into Stack
    int i ;
    for (i=0 ; i<10 ; i++) {
        cout << "pushing value " << i*i << " in stack" << endl ;
        s.push(i*i) ;
    } 
    // Count doubles in fifo
    cout << s.nitems() << " values in stack" << endl ; 

    Stack sclone(s) ; // clone s with new variable sclone. "sclone(s)" calls the default copy constructor. Just as when you create "Stack s" you call "Stack()", now you do "Stack sclone(s)" and call the default Stack(const Stack&). 
    s.inspect() ; // inspect the stack when it is full
    sclone.inspect() ; // inspect the copied sclone. This looks okay 

    // Read doubles back from fifo 
    while (!s.empty()) { 
        double val = s.pop() ; 
        cout << "popping value " << val << " from stack" << endl ; 
    } 

    s.inspect() ; // insepct stack when it is empty 
    sclone.inspect() ; // inspection of sclone still looks good 

    for (i=0 ; i<5 ; i++) {
        cout << "pushing value " << i*100 << " in stack" << endl ;
        s.push(i*100) ;
    } 

    s.inspect() ; // insepct stack when it is empty 
    sclone.inspect() ; // 3.3e) --> inspection of sclone shows that the first 5 elements of sclone are changed by the values that were pushed to s. It seems to be that the sclone copies all memory of s. Then if s pops out numbers it essentially only decreases its size that it points to. This way sclone still points to the whole stack. However if s starts to write new numbers again to its stack then it overwrites also sclone. 

    return 0 ;
}

// 3.3c)    They seem to have identical content after copying. 
// 3.3d)    After popping values from s, it seems that sclone is still the same after 
//          inspectation. 
// 3.3e)    Inspection of sclone shows that the first 5 elements of sclone are changed 
//          by the values that were pushed to s. It seems to be that the sclone copies 
//          all memory of s. Then if s pops out numbers it essentially only decreases 
//          its size that it points to. This way sclone still points to the whole stack. 
//          However if s starts to write new numbers again to its stack then it overwrites 
//          also sclone. 
// 3.3h)    It behaves all correctly. 
// 
// To run the files:
//      g++ -c stack.cc
//      g++ -c main.cc
//      g++ -o ouput main.o stack.o 
//      ./output
// 
// To switch between git versions: 
//      git checkout t_ex31 (you cannot add and commit this, but you first need to go to a branch, it is sort of a draft version that you then need to add to a branch) 
//      git checkout -b b_ex31 (for new branch you need the "-b" term with it) 
//      (edit, git add, git commit) (important to commit before switching branches!)
//      git checkout master 
//      git merge b_ex31 
//      git push origin b_ex31 (to push to gitlab) 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 

