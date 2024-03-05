#include <iostream>
#include <list>
#include <vector>
using namespace std ; 

int main() {
    int n = 100000 ; 
    vector<int> array ; //create list/vector
    vector<int>::iterator iter ; // create iterator
    iter = array.begin() ; // set iterator to begin of list/vector
    for (int i=0; i<n; i++) { 
        array.push_back(i) ; // put a number at end of array
        ++iter ; // go further with the iterator
        // cout << *iter << " " ; // if you want to print the value 
        }
    // cout << endl << endl ; // 

    for (iter=array.begin(); iter!=array.end(); ++iter) {
        if (*iter%3==0) {iter = array.erase(iter) ;}
    }

    // for (iter=array.begin(); iter!=array.end(); ++iter) { // to look if the third element is removed
    //     cout << *iter << " " ; // show value
    // }
    // cout << endl ; 

    return 0 ; 
}

// 7.3d)    With n=100.000 we got 
//          for "list":
//          0.01user 0.00system 0:00.01elapsed 100%CPU (0avgtext+0avgdata 6576maxresident)k
//          for "vector":
//          0.11user 0.00system 0:00.11elapsed 100%CPU (0avgtext+0avgdata 3556maxresident)k
//          So the list is about 10 times faster here 

