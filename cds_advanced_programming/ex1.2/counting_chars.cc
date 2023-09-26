#include <iostream>

int main() {
    char word[256] ; // make room for 265 characters
    char *word_pointer = &word[0] ; // create a pointer to the first character of the word-array


    std::cout << "please type something: " ;
    std::cin.getline(word, 256) ; // set word to the value that is typed in 
    std::cout << "memory address of first letter: " << &word_pointer << "\n" ; // get memory address of the word pointer

    std::cout << "print word with spaces: " ;
    while (*word_pointer != '\0') { //do while loop until end of word[256] which is closed by '\0'
        std::cout << *word_pointer << " " ; // print character where word pointer points to                 
        word_pointer++ ; // get the pointer to next character 
        }
    std::cout << "\n" ;


    word_pointer = &word[0] ; // set word pointer back to start 
    int upper_case = 0 ; // set upper case counter to 0
    int lower_case = 0 ; // set lower case counter to 0
    int spaces = 0 ; // set lower case counter to 0
    int special_char = 0 ; // set lower case counter to 0
    while (*word_pointer != '\0') { // loop over all character until end of the word
        if (int(*word_pointer) >= int('A') and int(*word_pointer) <= int('Z')) { // look at if integer value of the character is in the range of integers for the characters between upper case A and Z 
            upper_case += 1 ; 
        } else if (int(*word_pointer) >= int('a') and int(*word_pointer) <= int('z')) { // look at if integer value of the character is in the range of integers for the characters between lower case a and z
            lower_case += 1 ; 
        } else if (int(*word_pointer) == int(' ')) { // look at spaces
            spaces += 1 ; 
        } else if (int(*word_pointer) >= int('!') and int(*word_pointer) <= int('/')) { // look at special charachter segment in ASCII
            special_char += 1 ; 
        } else if (int(*word_pointer) >= int(':') and int(*word_pointer) <= int('@')) { // look at special charachter segment in ASCII
            special_char += 1 ; 
        } else if (int(*word_pointer) >= int('[') and int(*word_pointer) <= int('`')) { // look at special charachter segment in ASCII
            special_char += 1 ; 
        } else if (int(*word_pointer) >= int('{') and int(*word_pointer) <= int('~')) { // look at special charachter segment in ASCII
            special_char += 1 ; 
        } 
        word_pointer++ ; // get pointer to the next character
    }
    std::cout << "uppercase:          " << upper_case << "\n" ; // print amount of upper case letters counted
    std::cout << "lowercase:          " << lower_case << "\n" ;// print amount of lower case letters counted
    std::cout << "spaces:             " << spaces << "\n" ;// print amount of lower case letters counted
    std::cout << "special characters: " << special_char << "\n" ;// print amount of lower case letters counted


    return 0 ; 
}

