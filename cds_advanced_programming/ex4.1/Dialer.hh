#include <iostream>
#ifndef DIALER_HH
#define DIALER_HH
#include "Button.hh"

class Dialer {
public:

  Dialer() { 
    std::cout << "Dialer Constructor " << this << std::endl ; 
    buttons = new Button[12] ; // creates new space for the pointer 
  }
  
  // Dialer(const Dialer& dialer) { 
  //   std::cout << "Dialer Copy Constructor " << this << std::endl ; 
  //   buttons = new Button[12] ; // creates new space for the pointer 
  //   for (int i=0; i<12; i++) {
  //     *(buttons+i) = *(dialer.buttons+i) ; // sets the values of the array/pointer to the one you copy from 
  //   }
  
  Dialer(const Dialer& dialer) { 
    std::cout << "Dialer Copy Constructor " << this << std::endl ; 
    buttons = new Button[12] ; // creates new space for the pointer 
    for (int i=0; i<12; i++) {
      buttons[i] = dialer.buttons[i] ; // sets the values of the array/pointer to the one you copy from 
    }
  }

  ~Dialer() { 
    std::cout << "Dialer Destructor " << this << std::endl ; 
    delete[] buttons; // delete the allocated memory
  }

private:

  Button* buttons ; // create button pointer where you will store memory later 

} ;

#endif 

