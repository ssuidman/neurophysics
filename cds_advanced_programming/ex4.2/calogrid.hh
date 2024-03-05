#ifndef CALOGRID 
#define CALOGRID 

#include <iostream> 
#include "calocell.hh" 

class Calogrid { 
    public: 
        Calogrid(int nx0, int ny0) ; 
        Calogrid(const Calogrid& calogrid) ; 
        ~Calogrid() ; 
        Calocell* cell(int x, int y) ; 
        const Calocell* cell (int x, int y) const ; // this cannot use set_energy() for example, because it cannot change the energy 

    private: 
        int nx ; 
        int ny ; 
        Calocell** elements ; // pointer to array of pointers
        // Calocell elements[nx*ny] // This would give the error of 4.2n
} ;

#endif

