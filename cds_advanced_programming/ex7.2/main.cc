#include <iostream>
#include <vector>
#include <string>
using namespace std ; 

template<class T>
bool isPalindrome(vector<T> v) ;

int main() {
    vector<int> vec = {1,2,3,4,5,4,3,2,1} ;
    vector<int>::iterator iter1 ;
    vector<int>::iterator iter2 ;
    iter1 = vec.begin() ; // set iter1 to beginning
    iter2 = vec.end() ; // set iter2 to end 
    iter2-- ; // set iter2 to begin of last element
    bool palindrome = true ; 

    while(iter2!=vec.begin()) {
        if(*iter1 != *iter2) { 
            palindrome = false ;
            }
        ++iter1 ;
        --iter2 ;
        
    }
    cout << palindrome << endl ; 

    
    vector<string> vec2 = {"apple","butter","cream","butter","appl"} ;
    bool check = isPalindrome(vec2) ;
    cout << (check ? "This is a palindrome" : "This is not a palindrome") << endl ; 
    return 0 ; 
}

template<class T>
bool isPalindrome(vector<T> v) {
    typename vector<T>::iterator iter1 ;
    typename vector<T>::iterator iter2 ;
    iter1 = v.begin() ; // set iter1 to beginning
    iter2 = v.end() ; // set iter2 to end 
    iter2-- ; // set iter2 to begin of last element
    bool palindrome = true ; 

    while(iter2!=v.begin()) {
        if(*iter1 != *iter2) { 
            palindrome = false ;
            }
        ++iter1 ;
        --iter2 ;
    }
    return palindrome ;
}  