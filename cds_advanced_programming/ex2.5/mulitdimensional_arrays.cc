#include <iostream>

double* multiply(double ivec[], int N, double mtx[][3]) { 
    // NOTE: Caller takes ownership of the pointer "return_array"
    double* return_array = new double[3] ; // create new pointer (=array)
    for (int i=0; i<3; i++) {
        int counter = 0 ; // counter is to keep track of the sum
        for (int j=0; j<N; j++) {
            counter += mtx[j][i]*ivec[j] ; // add each time the elements of matrix multiplication
        }
        *(return_array+i) = counter ; //set the value of the array to this multiplication
    }
    return return_array ; // return the array (=pointer)
}


int main() {
    int N = 5; 
    double vector[N] = {1,2,3,4,5} ; // create array (=pointer)
    double matrix[N][3] = {{1,3,5},{7,9,0},{2,4,6},{8,1,3},{5,7,9}} ; // create matrix, which is pointer to different pointers (that are pointing to doubles)
    double* return_vector = multiply(vector,N,matrix) ; // create pointer for the result and take owner ship of memory!
    
    std::cout << "matrix multiplication result:    [ " ;
    for (int i=0; i<3; i++) {
        std::cout << *(return_vector+i) << " " ; // print values of the array after multiplication
    }
    std::cout << "]" << std::endl ; 

    for (int i=0; i<N; i++) {
        for (int j=0; j<3; j++) {
            std::cout << *(*(matrix+i)+j) << " " ; // just as with ex1.2, "matrix" is a pointer that points to other pointers. Now if you do (matrix+i) you switch between pointers to groups of numbers (so for example the pointer to {1,3,5}). Then to get the value of the first pointer to this group you do *(matrix+i) and to get the j'th pointer (for example the second(=j'th) pointer to the number "3" inside {1,3,5}) you do (*(matrix+i)+j). Then if you want to get the value of this pointer (which is 3 in this example) you do *(*(matrix+i)+j)
            std::cout << matrix[i][j] << " " ; // print matrix[i][j] the easy way instead of pointerlike
        }
    }
    std::cout << std::endl ;

    delete[] return_vector ; // delete memory that you have become owner of 
    return 0 ;
}

