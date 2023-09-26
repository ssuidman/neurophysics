#include <iostream>
#include "calorimeter.hh"

Calogrid& Calorimeter::grid() {
    return calogrid ;
}

const Calogrid& Calorimeter::grid() const {
    return calogrid ;
} 

Point& Calorimeter::position() {
    return point ;
}

const Point& Calorimeter::position() const {
    return point ;
}

