// Question:    Why does running this programming (the right way: "./output example.txt") give the
//              wrong amount of characters by a difference of exactly the amount of lines (80)?

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
using namespace std ;

int main(int argc, char* argv[] ) {
    ifstream ifs(*(argv+1)) ; // get file. argc is 2 and "*(argv+0)" is "./output", while "*(argv+1)" is the file name
    if(!ifs.is_open()) {cout << "File is not recognized" << endl ;} // check if file is opened

    int line_count = 0 ; // count lines
    int char_count = 0 ; // count char
    int word_count = 0 ; // count words
    char buf[100] ; // set the buffer for a line of 100 characters 

    while (ifs.getline(buf,100)) {  // if there are no more lines it is ended
        char_count += strlen(buf)+1 ; // 
        istringstream iss(buf) ; // set iss to buf at the moment 

        cout << strlen(buf) << " " ;
        while (iss >> buf) { // go from word to word in buf until you are at the end of the stream
            cout << buf << " " ; // printing the text
            word_count += 1;} // count the words
        line_count += 1 ; // count the lines

        cout << endl ; 
    } 
    cout << endl ;
    cout << line_count << " " << word_count << " " << char_count << endl ; // There are in total 80 lines (checked with "wc -l example.txt"). Also there are 4247 characters in total ("wc -c example.txt" gives exact 80 characters more, but online checking gives the same)
    
    return 0 ;
}

