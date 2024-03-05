#include <iostream>
using namespace std ;


char* memory_leak () {
    char* leaker = new char[10000] ; // each time there is a memory loss for 10000 charachters 
    delete leaker ; // when removing this line the %MEM goes up to 30 and then stops
    return 0 ;
}


int main() {
    int N = 10000000 ;
    for (int i =0; i<N; i++) {
        memory_leak() ; 
    }
    return 0 ;
}


// To look at the memory you can type "top -u ssuidman" in another terminal to look at the process. In the columns %MEM you can then see building up the memory loss. 
// When running with to much memory loss the shell connection is lost temporarily and starts again. The following error occurs:
// "Restarting the terminal because the connection to the shell process was lost...""