#include <iostream>
#include <list>
#include "Square.hh"
#include "Circle.hh"
using namespace std ; 

void listShapes(list<Shape*> l) ;

int main() {
    Square square(5) ; 
    cout << "Square surface = " << square.surface() << "; circumference = " << square.circumference() << endl ; 
    Circle circle(5) ; 
    cout << "Circle surface = " << circle.surface() << "; circumference = " << circle.circumference() << endl ; 

    list<Shape*> shapes ;
    Square square1(1) ;
    Square square2(2) ;
    Square square3(3) ;
    Circle circle1(1) ;
    Circle circle2(2) ;
    Circle circle3(3) ;
    shapes.push_back(&square1) ;
    shapes.push_back(&square2) ;
    shapes.push_back(&square3) ;
    shapes.push_back(&circle1) ;
    shapes.push_back(&circle2) ;
    shapes.push_back(&circle3) ;
    listShapes(shapes) ;

    return 0 ; 
}


void listShapes(list<Shape*> l) {
    list<Shape*>::iterator iter ; 
    for (iter=l.begin(); iter!=l.end(); ++iter) { 
        cout << (*iter)->shapeName() << ": surface = " << (*iter)->surface() << "; circumference = " << (*iter)->circumference() << endl ; 
    } 
} 

// 8.3f/g)  The code does not compile right, because of the following error: 
//          "
//          cannot declare variable ‘square’ to be of abstract type ‘Square’
//          because the following virtual functions are pure within ‘Square’:
//          ‘virtual const char* Shape::shapeName() const’ 
//          "
//          The problem is that shapeName() on itself cannot exist. There need 
//          to be shapes such as "circle" or "square" that also own this funcion. 
//          When adding these functions there the code compiles the right way. 
