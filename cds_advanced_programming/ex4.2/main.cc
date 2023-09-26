#include <iostream>
#include "calocell.hh"
#include "point.hh"
#include "calogrid.hh"
#include "calorimeter.hh"

int main() {
    Calocell calocell(0.1, 0) ; // initialize
    std::cout << calocell.get_energy() << std::endl ; 
    calocell.set_energy(1.1) ; // change value
    std::cout << calocell.get_energy() << std::endl ; 
    std::cout << std::endl ; 
    
    Point point ; // default (0,0,0)
    Point point2(1,2,3) ; // set own initial values 
    std::cout << "Point 1: (" << point.get_x() << "," << point.get_y() << "," << point.get_z() << ")" << std::endl ;
    std::cout << "Point 2: (" << point2.get_x() << "," << point2.get_y() << "," << point2.get_z() << ")" << std::endl ;
    point.set_xyz(4,5,6) ;
    point2.set_x(0) ;
    std::cout << "Point 1: (" << point.get_x() << "," << point.get_y() << "," << point.get_z() << ")" << std::endl ;
    std::cout << "Point 2: (" << point2.get_x() << "," << point2.get_y() << "," << point2.get_z() << ")" << std::endl ;
    std::cout << std::endl ; 

    Calogrid calogrid(2,3) ; // call the constructor 
    Calogrid calogrid2(calogrid) ; // call the copy constructor
    Calocell* place = calogrid.cell(1,2) ; // NOT "const Calocell* place", because then set_energy cannot be used
    place->set_energy(2.1) ;  // set energy of cell
    std::cout << place->get_energy() << std::endl ; // use the "->"-operator because place is a pointer! This does the same as *(place).set_energy(2.1)
    std::cout << std::endl ; 
    
    Calorimeter calorimeter(1,3) ; // create calorimeter instance with standard coordinates (0,0,0)
    Calorimeter calorimeter2(2,1,1.1,2.5,-0.2) ; // create calorimeter of 2x1 with self set coordinates 
    Point& point3 = calorimeter2.position() ; // get point of calorimeter 
    const Point& point4 = calorimeter2.position() ; //just to let know that this can also be constructed, but here you cannot use get_x(), etc
    Calogrid& calogrid3 = calorimeter2.grid() ; // get calogrid (reference)
    Calocell* place2 = calogrid3.cell(1,0) ; // get (pointer to) certain place cell (within the calogrid, otherwise segmentation fault)
    Calocell* place3 = calorimeter2.grid().cell(0,0) ; // call a place cell at once from the calorimeter 
    place2->set_energy(3.3) ; // set energy of the place cell to 3.3
    place3->set_energy(420) ; // set energy of other place cell
    point3.set_x(12.1) ; // set x value of created point 
    std::cout << "Point 3: (" << point3.get_x() << "," << point3.get_y() << "," << point3.get_z() << ")" << std::endl ; //print coordinates of the place cell 
    std::cout << place2->get_energy() << " " << place3->get_energy() << std::endl ; // print the energy of this cell that you just set. 

    return 0 ;


}

// 4.2b)    Copy constructor not needed, because you won't need to copy the cell. 
//          Destructor also not, because there is no dynamic memory allocation that 
//          needs to be deleted. 
// 
// 4.2d)    Const is used with the get_energy() and get_ID(), because you only want to
//          call energy and ID, but keep them unchanged. 
// 
// 4.2g)    You don't need a destructor, because there is no dynamic memory allocation, 
//          but you need a copy constructor, because if you want to change a variable 
//          point_A, but also want to know the old coordinates you need a point_B that 
//          copies point_A. You can of course argue that you can create a new point_B,
//          instead of changing point_A, but it does not hurt to create a copy 
//          constructor anyway. 
// 
// 4.2u)    Copy constructor and destructor are not necessary, because you don't want
//          to make another calorimeter and are not allocating new memory. 