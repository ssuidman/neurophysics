#include <iostream>
#include <cstring>
using namespace std ;



void order_1(int& a, int& b) { 
    int x ; // create storing integer
    if (a > b) { // when you replace ">" by "<" the algorithm sorts from bit to small instead of small to big
        x = a ; // temporary int x is set to the a outside of the function
        a = b ; // the a outside of the function is replaced 
        b = x ; // the b outside of the function is replaced 
    }
}

int* sort_1(int* array) { // function that sorts numbers 
    // NB: Caller take ownership of *sorted !
    int n = 10 ; // number of integers to be sorted
    int* sorted = new int[n] ; // create new pointer to return out of this function, in sort_3 this is done without creating new pointers

    for (int i=0; i<n; i++) { 
        *(sorted+i) = array[i] ; // set all values of sorted (return pointer) to the array values that you have from input without shifting the pointer
    } 

    // sorted = zero_ptr ; // set sorted pointer back to start
    for (int i=0; i<n-1; i++) { // use algorithm that has been given
        for (int j=0; j<n-1-i; j++) { 
            order_1(*(sorted+j),*(sorted+j+1)) ; // swap the two pointers to the values by giving the values of the pointers themselves
        } 
    }

    return sorted ;  // return pointer out of function
} 



void order_2(int* a, int* b) {
    // Caller takes ownership of swap! 
    int x ; 
    if (*a > *b) { // when you replace ">" by "<" the algorithm sorts from bit to small instead of small to big
        x = *a ; // setting x to the value of the pointer a
        *a = *b ; // setting the value of pointer a outside the function to the value of pointer b outside the function
        *b = x ;  // setting the value of pointer b outside the function to the value of the temporary pointer x
    }
}

int* sort_2(int* array) { // function that sorts numbers 
    // NB: Caller take ownership of *sorted !
    int n = 10 ; // number of integers to be sorted
    int* sorted = new int[n] ; // create new pointer to return out of this function

    for (int i=0; i<n; i++) { 
        *(sorted+i) = array[i] ; // set all values of sorted (return pointer) to the array values that you have from input without shifting the pointer
    } 

    for (int i=0; i<n-1; i++) { // use algorithm that has been given
        for (int j=0; j<n-1-i; j++) { 
            order_2((sorted+j),(sorted+j+1)) ; // swap the values that the pointers are referring to 
        } 
    }
    return sorted ;
}



void main_2() ;

int main() {
    int random[10] = {534,29,7,5111,703,119,90,1502,98,36} ; // create random numbers yourself
    int* sorted_array_1 = sort_1(random) ; // need to delete sorted_output afterwards
    int* sorted_array_2 = sort_2(random) ; // need to delete sorted_output afterwards

    for (int i=0; i<10; i++) {
        cout << *(sorted_array_1+i) << " " ; // print the values inside the pointer
        }
    cout << endl ;

    for (int i=0; i<10; i++) {
        cout << *(sorted_array_2+i) << " " ; // print the values inside the pointer
        }
    cout << endl ;

    delete[] sorted_array_1 ; // delete sorted_output pointer so there is no memory leakage 
    delete[] sorted_array_2 ; // delete sorted_output pointer so there is no memory leakage 

    main_2() ; // run the adapted main function

    return 0 ; 
} 



void order_3(const char** a, const char** b) { // takes two pointers to an array of pointers 
    const char* x ; // create temporary pointer 
    int comparer = strcmp(*a,*b) ; // compares two pointers to (full) strings and returns "0","<0",",">0"-numbers
    
    if (comparer>0) { // when you replace ">" by "<", the algorithm sorts from small to big instead of big to small 
        x = *a ; // set temporary pointer to pointer to string a (that is outside of the function)
        *a = *b ; // set pointer to string b (that is outside of the function) to pointer to string a (that is outside of the function)
        *b = x ; // set pointer to string b (that is outside of the function) to temproary pointer
    }

}

void sort_3(const char** array, int n) { // the char** is a pointer that points to an array of pointers
    for (int i=0; i<n-1; i++) { // use algorithm that has been given
        for (int j=0; j<n-1-i; j++) { 
            order_3((array+j),(array+j+1)) ; // swap the order that the pointer points to the constant string-pointers
        } 
    }
}



void main_2() { 
    int n = 10 ; // set the value of the amount of strings that you want to create
    const char* string_array[n] = {"appel","kaas","A","Z","cholera","youtube","++c","touw","regels","nalaten"} ; // create a pointer to an array of pointers to strings. To get the pointer to the second (full) string you can do *(string_array+1) for example. The pointers to the strings are constant, but the pointer to the pointers is not!
    sort_3(string_array,n) ; // the pointer that points to the other pointers is taken as input (as is n=10) and the order the pointer points will be swapped

    for (int i=0; i<10; i++) {
        cout << *(string_array+i) << " " ; // print the sorted string_array
        }
    cout << endl ;
} 

// 2.1f)    You need the const here, because you want to order the words, but not the letters 
//          in the words. So you want to keep the words themselves constant. 
// 
