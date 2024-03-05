#include <iostream>

template<class T>
T* sort(T* array, int n) ; // sorting algorithm

template<class T>
void order(T& a, T& b) ; // ordering algorithm

template<>
void order(const char*& a, const char*& b) ; // ordering algorithm for const char*

template<class T>
void display(T t[], int n) ; // display algorithm



template<class T>
T* sort(T* array, int n) { // NB: Caller take ownership of *sorted !
    T* sorted = new T[n] ; // create new pointer to return out of this function, in sort_3 this is done without creating new pointers

    for (int i=0; i<n; i++) { // set values of pointer to array values 
        *(sorted+i) = array[i] ; 
        } 

    for (int i=0; i<n-1; i++) { // use algorithm that has been given
        for (int j=0; j<n-1-i; j++) { 
            order(*(sorted+j),*(sorted+j+1)) ; // swap the two pointers to the values by giving the values of the pointers themselves
        } 
    }

    return sorted ;  // return pointer out of function
} 

template<class T>
void order(T& a, T& b) { 
    T x ; // create storing integer
    if (a > b) { // when you replace ">" by "<" the algorithm sorts from bit to small instead of small to big
        x = a ; // temporary int x is set to the a outside of the function
        a = b ; // the a outside of the function is replaced 
        b = x ; // the b outside of the function is replaced 
    }
}

template<>
void order(const char*& a, const char*& b) { // takes two references to (word) pointers
    const char* x ; // create a constant char pointer 
    if (strcmp(a,b)>0) { // compare two pointers to words 
        x = a ; // swapping
        a = b ;
        b = x ;
    }
}

template<class T>
void display(T t[], int n) { // t is a pointer to class T
    for (int i=0; i<n; i++) {
        std::cout << *(t+i) << " " ; // print the values inside the pointer
        }
        std::cout << std::endl ; 
}


