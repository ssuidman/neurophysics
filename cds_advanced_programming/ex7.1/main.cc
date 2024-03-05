#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
using namespace std ; 

int main() {
    ifstream ifs("../ex5.3/example.txt") ; 
    string word ; // make string for the words
    map<string,int> MyMap ; // make a map to store the strings
    map<string,int>::iterator iter ; // make an iterator to go through the map

    while(getline(ifs,word)) { // printing the example.txt while putting the words in My Map with counts. 
        istringstream line(word) ;
        while(line >> word) { 
            cout << word << " " ; 
            MyMap[word] += 1 ; // put the words in the map (if they are not yet in it) and raise the word counter by 1 (with zero as default parameter). 
            }
        cout << endl ; 
    }
    // ifs is now at the end of the text, so if you would start over again you need to set ifs to the beginning of the text. 

    iter = MyMap.begin() ; // set iterator to begin of the map
    while(iter != MyMap.end()) { // if iterator at the end --> stop
        cout << iter->first << " " << iter->second << endl ; // print the words with the counts
        ++iter; // go to the next "MyMap" element
    }

    return 0 ; 
}


