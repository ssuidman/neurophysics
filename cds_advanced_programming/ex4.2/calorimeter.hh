#ifndef CALORIMETER
#define CALORIMETER

#include <iostream>
#include "point.hh"
#include "calogrid.hh"

class Calorimeter {
    public:
        Calorimeter(int nx0, int ny0) : nx(nx0), ny(ny0), x(0), y(0), z(0), calogrid(nx0, ny0), point(0,0,0) {std::cout << "Calorimeter constructor" << std::endl ;} // constructor of calorimeter with inital values set
        Calorimeter(int nx0, int ny0, double x0, double y0, double z0) : nx(nx0), ny(ny0), x(x0), y(y0), z(z0), calogrid(nx0, ny0), point(x0,y0,z0) {std::cout << "Calorimeter constructor" << std::endl ;}
        Calogrid& grid() ; 
        const Calogrid& grid() const ; 
        Point& position() ; 
        const Point& position() const ; 
    private: 
        int nx ;
        int ny ;
        double x ;
        double y ;
        double z ;
        Calogrid calogrid ;
        Point point ; 
} ;

#endif

