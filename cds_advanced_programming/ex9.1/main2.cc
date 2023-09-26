#include <iostream>
#include <thread>
#include <string>
#include <chrono>
using namespace std;

void f1() { 
    cout << "f1 first " << endl ;
    this_thread::sleep_for(chrono::seconds(2));
    cout << "f1 second " << endl ;
}
  
void f2(const string& s) { 
    this_thread::sleep_for(chrono::seconds(1));
    cout << s << endl ;
}

void f3(int i) {
    cout << i << endl ; 
}

void f4(int& i) {
    cout << i << endl ; 
}

int main() {
    thread t1(f1);
    thread t2{f2, "f2 first "};
    t1.join(); // deleting this line gives a run-time error ("terminate called without an active exception" and "Aborted")
    t2.join();

    int a = 1 ;
    int& b = a ;
    thread t3(f3,a) ; // passing by int
    // thread t4(f4,b) ; // this does not work
    thread t4(f4,ref(a)) ; // passing by reference using "std::ref(a)"
    thread t5(f3,10) ; // passing by value 
    t3.join();
    t4.join();
    t5.join();
    return 0;
}

// General:     Run the program with "g++ -o output main.cc -pthread", otherwise errors. 
// 9.1d)        You need to use std::ref() to get a reference passed in a thread. Passing 
//              by value is also easy. 

