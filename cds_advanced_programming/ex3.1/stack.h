
#ifndef STACK
#define STACK

const int LEN = 80 ; // default stack length

class Stack {     // a class is by default private while a struct is by default public 

  public:   // Interface (=public normally)
    // void init() { count = 0 ; } // not public any more, because of constructor
    Stack() { init() ; } // constructor --> this self initializes the function
    ~Stack() { close() ; } // destructor --> removes variables that were initialized
    int nitems() { return count ; }
    bool full() { return (count==LEN) ; }
    bool empty() { return (count==0) ; }

    void push(double c) ;
    double pop() ;
    void inspect() ;

  private:  // Implementation (=private normally)
    void init() { count = 0 ; std::cout << "initialized by constructor" << std::endl ;} // after the constructor was added, init() is now private, so there can only be self initialization
    void close() { std::cout << "closed by destructor" << std::endl ; } // how to remove count and when is removed
    double s[LEN] ;
    int count ; // --> why can you put the definition of count after when it is used? 
} ;

#endif


