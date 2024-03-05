#include <iostream>

namespace Black {  // namespace Black is created
  void print(int k) {} ; 
} 

namespace White { 
  void print(int k) {} ; 
} 

// using Black::print ; // b) Global using declaration -- OK? -----> not okay to use, because in "void print(int k)" the function print(k-1 is called. This kan either be itself (so "void print(int k)" or "Black::print" (if you write this line). 
 
void sub1() { 
  using White::print ; // Local using declaration 
  std::cout << "sub1: " << 5 << std::endl ; // look if this function is used by showing the number 5 through the usage of the std namespace from the standard library
  print(5) ;           // a) Which print() is called? ---> the print from the namespace White is used 
} 

void print(int k) { 
  if (k>0) { 
    std::cout << "print: " << k << std::endl ; // look at what functions are used 
    print(k-1) ;       // a) Which print() is called? --> This function itself ("void print(int k)") is used, because there is no explicit declaration of a namespace. 
  } 
} 

int main() {
    sub1() ; // call the functions sub1 and print to see what happens
    print(10) ;
    return 0 ;
}