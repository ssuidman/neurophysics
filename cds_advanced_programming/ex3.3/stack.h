#ifndef STACK
#define STACK


class Stack {     // a class is by default private while a struct is by default public 

  public:   // Interface (=public normally)
    Stack(int length) { init(length) ;} // Take ownership of "s". It can be good to specify length, because if you want a big stack and know this already you don't need to call the grow function all the time. 
    Stack() { init() ;} // constructor --> take ownership of "s".
    ~Stack() { close() ; } // destructor --> removes "s" that was initialized in the beginning 
    Stack(const Stack& s_old) ; // The input is const and then Stack &s_olc. This is just as you would have a function f(int i) only now the input is of the class you created yourself "Stack" 
    int nitems() { return count ; } 
    bool full() { return (count==LEN) ; }
    bool empty() { return (count==0) ; }
    void push(double c) ;
    double pop() ;
    void inspect() ;
    void grow(int delta) ;
    double* sclone() ;

  private:  // Implementation (=private normally)
    void init(int length) ;    
    void init() ;    
    void close() ;
    // void close(const Stack& s) ;
    int LEN ; // overall LEN variable that is used and is initialized in the init() function
    double* s ;
    int count ; // --> why can you put the definition of count after when it is used? 
} ; 

#endif

