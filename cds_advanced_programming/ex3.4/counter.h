#ifndef COUNTER 
#define COUNTER 

class Counter {
    public:
        Counter() { count++ ; } // if you create (constructor) new instance the static int count goes up by 1
        ~Counter() { count-- ; } // if you destruct instance the static int count goes down by 1
        static int getCounter() ; // declare the static function that returns the counter

    private:
        static int count ; // declare the static int count variable 
} ;

#endif 

