#include <iostream>
#include "calogrid.hh"

Calogrid::Calogrid(int nx0, int ny0) {
    std::cout << "Calogrid constructor " << this << std::endl;
    nx=nx0; 
    ny=ny0;
    elements = new Calocell*[nx*ny]; // create new pointer to pointers with length nx*ny
    double energy0 = 0.4 ;
    int ID0 = 3 ;
    for (int i=0; i<nx*ny; i++) {
        *(elements+i) = new Calocell(energy0, ID0) ;
    }
}

Calogrid::Calogrid(const Calogrid& calogrid) {
    std::cout << "Calogrid copy constructor " << this << std::endl;
    nx=calogrid.nx; 
    ny=calogrid.ny;
    elements = new Calocell*[nx*ny];
    for (int i=0; i<nx*ny; i++) {
        *(elements+i) = *(calogrid.elements+i) ;
        }
}

Calogrid::~Calogrid() {
    std::cout << "Calogrid destructor " << this << std::endl;
    delete[] elements ;
}

Calocell* Calogrid::cell(int x, int y) {
    Calocell* certain_cell ;
    if (x>=0 and x<nx and y>=0 and y<ny) {
        certain_cell = elements[x+y*nx]; 
    } else {
        certain_cell = NULL ;
    }
    return certain_cell ;
}

const Calocell* Calogrid::cell (int x, int y) const {
    const Calocell* certain_cell ;
    if (x>=0 and x<nx and y>=0 and y<ny) {
        certain_cell = elements[x+y*nx]; 
    } else {
        certain_cell = NULL ;
    }
    return certain_cell ;
}

