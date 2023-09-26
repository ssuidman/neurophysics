#ifndef SQUARE_HH
#define SQUARE_HH

#include "Shape.hh"

class Square: public Shape {
public:

  // Constructor, destructor
  Square(double size) : _size(size) {} ;
  virtual ~Square() {} ;

  // Implementation of abstract interface
  virtual double surface() const { return _size * _size ; }
  virtual double circumference() const { return 4 * _size ; }
  virtual const char* shapeName() const { return "Square" ; }

private:

  double _size ;

} ;

#endif
