#include <iostream>
#include "Array.hh"
using namespace std ; 

int main() {
    Array<int> integers(10,7) ;

    const char default_val = 'a' ; // default const char, where you need to take reference from
    const char set_val = 'b' ; // one value
    Array<const char*> characters(8,&default_val) ; // set default const char*

    integers[1000] = 12 ; // set a value larger than size of "integers"
    Array<int> reproduce(1,1); 
    reproduce = integers ; // use the "=" operator
    
    for (int i=998; i<1003; i++) {
        cout << integers[i] << " " ; // print out some integers
    }
    cout << endl ; 

    for (int i=998; i<1003; i++) {
        cout << reproduce[i] << " " ; // print out the reproduced array. You can see it going wrong after 1000, because there the array produces 1's (it's initial values). 
    }
    cout << endl ; 

    characters[1000] = &set_val ; // set value of const char*
    for (int i=998; i<1003; i++) {
        cout << *characters[i] << " " ; // print out the value of "const char*" using the "*"-symbol
    }
    cout << endl ; 



    return 0 ;
}

// 6.2d)    If you try to get the 1000'th element of a 10 sized array, you 
//          can still change the value. However the size of the array is not 
//          changed, so it does not look safe. 
