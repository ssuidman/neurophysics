#include <iostream>
#include "counter.h"
using namespace std ;

int Counter::count = 0 ; // set the static variable to zero outside the declaration. 

int Counter::getCounter() { // return the counter 
    return count ;
}

