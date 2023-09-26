#ifndef STACK
#define STACK


class Stack {     // a class is by default private while a struct is by default public 

  public:   // Interface (=public normally)
    Stack(int length) { init(length) ;} // Take ownership of "s". It can be good to specify length, because if you want a big stack and know this already you don't need to call the grow function all the time. 
    Stack() { init() ;} // Take ownership of "s".
    ~Stack() { close() ; } // destructor --> removes "s" that was initialized in the beginning 
    int nitems() { return count ; }
    bool full() { return (count==LEN) ; }
    bool empty() { return (count==0) ; }
    void push(double c) ;
    double pop() ;
    void inspect() ;
    void grow(int delta) ;

  private:  // Implementation (=private normally)
    void init(int length) ;    
    void init() ;    
    void close() ;
    int LEN ; // overall LEN variable that is used and is initialized in the init() function
    double* s ;
    int count ; // --> why can you put the definition of count after when it is used? 
} ; 

#endif

