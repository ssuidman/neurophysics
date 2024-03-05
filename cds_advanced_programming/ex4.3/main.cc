#include <iostream>
#include "String.hh"
#include <cstring>
using namespace std ;


String operator+(const String& s1, const String& s2) { 
    String result(s1) ; // use copy constructor to make a total string that contains the added s1 (and in the end also s2)
    result += s2 ; // add s2 to s1 with the "+=" operator that is already created
    return result ; 
} 


int main() {
    String string("abc"), string2("defg"), string3, string4("hijkl"), string5 ;
    string4 += string2 ;
    string5 = string + string4 ;
    
    String s("Blah") ;   
    s += "Blah" ;

    string3 = string2 = string ; 
    cout << string.data() << " " << string2.data() << " " << string3.data() << " " << string4.data() << " " << string5.data() << " " << s.data() << endl ;
    cout << strlen(string) << endl ;
    return 0 ;
}

// 4.3i)    It is possible to add "Blah" and "Blah", because the compiler is converting 
//          a (const char*) to a string, hence you add two strings. 