#include <iostream> 
#include <string.h> // get stringlen package
using namespace std ; 


char* join1(const char* char1, const char* char2) {
    // NOTE: caller takes ownership of memory
    int char1_len = strlen(char1) ;
    int char2_len = strlen(char2) ;
    char* join1_pointer = new char[char1_len+char2_len+1] ; // creates pointer with a space of the total amount of characters 
    char* join1_pointer2 = join1_pointer ; // set  another pointer to the beginning of this space
    for (int i=0; i<char1_len; i++) { // loop over first word
        *join1_pointer = char1[i] ; // appoint word of character to new space where pointer is pointed to
        join1_pointer++ ; // get pointer to next (empty) character 
        // cout << *joint_pointer ; 
    } 
    for (int i=0; i<char2_len; i++) {
        *join1_pointer = char2[i] ; // same for second word
        join1_pointer++ ;
        // cout << *joint_pointer ; 
    } 
    *join1_pointer = '\0' ;
    join1_pointer = join1_pointer2 ;
    return join1_pointer ;
}


char* join2(const char* char1,const char* char2) { // same but easier function as the precious one by using strcat
    // NOTE: caller takes ownership of memory
    int char1_len = strlen(char1) ;
    int char2_len = strlen(char2) ;
    char* join2_pointer = new char[char1_len+char2_len+1] ;
    *join2_pointer = '\0' ;
    strcat(join2_pointer,char1) ;
    strcat(join2_pointer,char2) ;
    return join2_pointer ;
}

char* joinb(const char* char1,const char* char2) { // now using a space in between words
    // NOTE: caller takes ownership of memory
    int char1_len = strlen(char1) ;
    int char2_len = strlen(char2) ;
    char* joinb_pointer = new char[char1_len+char2_len+2] ;
    *joinb_pointer = '\0' ;
    strcat(joinb_pointer,char1) ;
    strcat(joinb_pointer," ") ;
    strcat(joinb_pointer,char2) ;
    return joinb_pointer ;
}

int main() {   
    
    char* join1_ptr = join1("alpha","bet") ; // print everything
    char* join2_ptr = join2("alpha","bet") ;
    char* joinb_ptr = joinb("duck","soup") ;

    cout << join1_ptr << endl ;   
    cout << join2_ptr << endl ;   
    cout << joinb_ptr << endl ;  
    
    delete[] join1_ptr ; // delete all pointers of used functions that were create by the new operator 
    delete[] join2_ptr ;   
    delete[] joinb_ptr ;

    return 0 ; 
} 

// Use "valgrind ./join_strings.out" to check for memory loss

// 1.3c)    Return type of join is a pointer to a character. 
// 1.3d)    You need to allocate new memory, because you want to store the full word in 
//          memory. To do that you need to allocate mememory for the characters. 
// 1.3e)    You need to allocate memory for the NULL character, hence the +1 in the code. 