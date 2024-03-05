#include <iostream>
#include <string.h>
using namespace std ;

int main() {
    unsigned int number ; // Can max fill in the number 2^32-1 = 4294967295
    cin >> number ; // get number from line 
    cout << "Your number is " << number << " and has a size of " << sizeof(number) << " bytes. " ; //sizeof(number) give the amount of bytes of a number 

    // Another way to do the following, is by using the ">>" operator as explained in the slides. 
    short int digit_holder[7]; // create array that holds binary number of 5 digits 
    short int compare = 0b11111 ;    // comparing number for five digits each time of input number. 5 is chosen bc 2^5
    short int k = 32 ; // 2^5 is division factor such that you "shift" each time 5 places in binary to the right 
    for (short int i = 0; i<7; i++) { // get five binary digits each time out of the number
        digit_holder[i] = number & compare ; // 01101 & 01011 = 01001 for example makes sure that the five digits of your number are represented in decimal and the rest in neglected of the number (vb after ... iterations is 1011011 & 11111 = 0011011 = 27)
        number = (number-digit_holder[i])/k ; // shift the binary number to the right by substracting terms that are already accounted for and then divide by 32. 
    } ;

    char base_32_table[32] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V'} ; // create 32-base table
    cout << "In 32-base it is represented by: " ;
    for (short int i = 0; i<7; i++) {
        int j = digit_holder[6-i] ; // print the 32-base number the right way from first to last and not otherwise. 
        cout << base_32_table[j] ;
    } ;
    cout << "." << endl ;
    return 0 ;
}

// Question:    The advantages of 16 base printing is for dealing with bit patterns in the 
//              computer. So 100111101111 can be dealth with more easily with hex or oct. 

