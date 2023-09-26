#include <iostream>
#include <iomanip>
using namespace std ;

int main() {
    int number ; 
    float f1, f2, f3 ;
    cin >> hex >> number ;
    cout << dec << number << endl << hex << number << endl << oct << number << endl ;

    cin >> f1 ;
    cin >> f2 ;
    cin >> f3 ;
    
    cout << scientific << f1 << " " << f2  << " " << f3 << endl ; 
    
    cout << setw(20) << "Value A     " << setw(20) << "Value B     " << setw(20) << "Value C     " << endl ;
    cout << setfill('-') << setw(61) << " " << setfill(' ') << endl ; 
    cout << setw(20) << f1 << setw(20) << f2  << setw(20) << f3 << endl ; 
    
    cout << fixed << left << setw(20) << setprecision(3) << f1 << setw(20) << setprecision(3) << f2 << setw(20) << setprecision(3) << f3 << endl ; 
    return 0 ;
}

// 5.1e)    Scientific does only need to be declared once. If you want to change this you 
//          can do did wherever you want. You can for example putting "fixed" right before 
//          "... << f3 ..."
// 
// 5.2f)    The term "setw(20)" needs to be before each number such that the number and the 
//          space are together 20. The text is printed on the right this way. 
// 
// 5.2g)    Putting "fixed" makes sure you switch from scientific to fixed mode. Now putting
//          "left" makes sure the number is printed on the left of your 20 digits field. 