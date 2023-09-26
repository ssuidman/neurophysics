#include <thread>
#include <string>
#include <iostream>
using namespace std;

void f1() { 
  cout << "first ";
  cout << endl ; // putting this line here results in printing the wrong order, such as "first third {new line} {new line} second" 
  cout << "second ";
}
  
void f2(const string& s) { 
  cout << s << endl ;
}

int main() {

  thread t1(f1);
  thread t2{f2, "third "};
  t1.join(); // deleting this line gives a run-time error ("terminate called without an active exception" and "Aborted")
  t2.join();
  
  return 0;
}

// General:   Run the program with "g++ -o output main.cc -pthread", otherwise errors. 
// 9.1a)      Putting the line "cout << endl ;" in f1() gives random order printing. 
// 9.1b)      deleting join gives run-time errors. 

