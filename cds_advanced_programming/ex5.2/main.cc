#include <iostream>
#include <fstream>
using namespace std ;

int main() {

    ifstream ifs("data1.txt") ; 
    cout << ifs.is_open() << endl ; // checks if file is open
    
    int number ; 
    while(ifs >> number) { // file ends if ifs is not set to a number any more 
        cout << number << endl ; 
        if (number==0) {break;} // to stop the printing when noticing a zero. 
    }
    return 0 ; 
}
