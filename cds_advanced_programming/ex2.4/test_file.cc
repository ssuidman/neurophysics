#include <iostream> // include build in header file with functions 
#include "min.hh" // include the min.hh file with all the function defined in it
#include "max.hh" // include the max.hh file with all the function defined in it

int min_func(int a,int b) { // int min function
    std::cout << "int function main file" << std::endl ;
    int minimum ; 
    if (a<=b) {
        minimum = a ;
    } else {
        minimum = b ;
    }
    return minimum ;
}

int main() {
    int a = 1 ;
    int b = 2 ;
    double c = mylib::min_func(a,b) ; // calling the int function in the library mylib in declared in min.hh and defined and in min.cc
    double c2 = min_func(a,b) ; //calling the int function in this file 
    int d[3] = {1,3,2} ;
    int e = mylib::max_func(d,3) ;
    std::cout << c << std::endl ;
    std::cout << c2 << std::endl ;
    std::cout << e << std::endl ;
    return 0 ;
}

// Two ways to couple everything together:
// 
// i)       Do after eachother "g++ -c min.cc", "g++ -c min.cc", "g++ -c library.cc". 
//          Then compile total file: "g++ -o total_file min.o max.o library.o" and run "./total_file"
// 
// ii)      With libraries you can also couple things together. 
// 
// 2.4d)    Look at "min.hh" and "max.hh" files. 
// 
// 2.5e)    Use libraries:
//          run "g++ -c min.cc" to create min.o to have compiled min file 
//          run "g++ -c max.cc" to create max.o to have complied max file
//          run "ar q libMinMax.a min.o max.o" to create library with min 
//          and max in it run "g++ -o test_file.o test_file.cc -L. libMinMax.a" 
//          to create complied test_file.o using the file test_file.cc and 
//          the library liMinMax.a where "L." give the pwd-folder to look 
//          for the library. If you want to update a file and the library, 
//          you run "g++ -c min.cc" then "ar r libMinMax.a min.o max.o", then
//          "g++ -o test_file.o test_file.cc -L. libMinMax.a" and then end
//          with "./test_file.o" to run the file. 
// 
// 2.5f)    There should be errors, because you have two functions that are 
//          almost exactly the same. I don't get this error, but it is still
//          important to define properly what function you want to call. 
// 
// 2.5g)    in all min.cc, max.cc, min.hh, max.hh files the function are put 
//          in the namespace mylib. 
// 
// 2.5h)    both int functions from the library and not the library are called. 
