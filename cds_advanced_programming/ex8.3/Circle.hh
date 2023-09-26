#ifndef CIRCLE_HH
#define CIRCLE_HH

#include "Shape.hh"
#include <cmath>

class Circle: public Shape {
public:

  // Constructor, destructor
  Circle(int radius) : _radius(radius) {} ;
  virtual ~Circle() {} ;

  // Implementation of abstract interface
  virtual double surface() const { return M_PI*_radius*_radius ; }
  virtual double circumference() const { return 2*M_PI*_radius ; }
  virtual const char* shapeName() const { return "Circle" ; }

private:

  int _radius ;

} ;

#endif
