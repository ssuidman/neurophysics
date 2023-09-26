#include <iostream>
#include <cstring>
#include "sort.hh"
using namespace std ;

int main() { 
    int rand_int[10] = {534,29,7,5111,703,119,90,1502,98,36} ; // create random numbers yourself 
    int* sort_int = sort(rand_int, 10) ; // need to delete sorted_output afterwards 
    float rand_float[8] = {4.3,8,30.2,9.1,-0.3,2.1,9.12,4814.2} ; 
    float* sort_float = sort(rand_float, 8) ; 
    const char* rand_char[9] = {"appel","kaas","A","Z","youtube","++c","touw","regels","nalaten"} ; // create a pointer pointing to word pointers. The pointers to words are const, because you don't want the words themselves changed
    const char** sort_char = sort(rand_char,9) ; // get the new pointer that is sorted and points to (const) word pointers

    display(sort_int, 10) ; 
    display(sort_float, 8) ; 
    display(sort_char, 9) ; 

    delete[] sort_int ; // delete sorted_output pointer so there is no memory leakage 
    delete[] sort_float ; 
    delete[] sort_char ; // you need to delete the memory that you have allocated for a new pointer that points to the right order of the word pionters. 

    return 0 ; 
} 

// 6.1h)    When using the order function it does not compares the characters 
//          themselves, but the values of the reference. So for example it 
//          compares 0x5640c2b32330 with 0x5640c2b32338. 
// 
// 6.1i)    The template function of sort can be used by different arrays, where:
//          T=int, T=foat or T=const char*. For "order" you need to make a specific 
//          function to compare the characters. 
