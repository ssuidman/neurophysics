#ifndef SHAPE_HH
#define SHAPE_HH

class Shape {
public:

  // Constructor, destructor
  Shape() {} ;
  virtual ~Shape() {} ;

  // Pure virtual interface functions
  virtual double surface() const = 0 ;
  virtual double circumference() const = 0 ;
  virtual const char* shapeName() const = 0 ;
} ;

#endif
