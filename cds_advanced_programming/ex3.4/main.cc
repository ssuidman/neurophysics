#include <iostream>
#include "counter.h"
using namespace std ;

int main() { 
  Counter a ; // create Counter instance
  Counter b ; // create another Counter instance
  cout << "there are now : " << Counter::getCounter()  << " Counter objects" << endl ; //Counter::get_Counter() returns count which is int
  if (true) { 
     Counter c ; // only inside the if statement this instance exists
     cout << "and now " << Counter::getCounter() << endl ; 
  } 
  cout << "and now " << Counter::getCounter() << endl ; // so now the instance "c" does not exist anymore 
} 

// 3.4b)  If you would initialize the int with the class declaration, each instance would initialize
//        this variable. You only want it to be initialized once, such that it can count all instances.
//        So it doesn't belong to an instance, but to the class itself. 