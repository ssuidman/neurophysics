#ifndef FUNCTIONS
#define FUNCTIONS
#include <iostream>
#include <string>
#include <queue>
#include <thread>
#include <chrono>
#include <mutex>
#include <condition_variable>
#include <random>
#include <pybind11/pybind11.h>
using namespace std ; 

class Threading { 
    public: 
        Threading (double time = 1.0) : t_max(1000*time) {cout << "constructor" << endl ; } 
        void run() ;

    private: 
        void produce() ; 
        void consume() ; 
        queue<double> q ; // create queue from class T (in this case double)
        bool condition ; 
        int t_max ; // in milliseconds
        mutex m ; 
        condition_variable cv ; 
        int i = 0 ; 
} ;


void Threading::run() { 
    thread t1 {&Threading::produce,this} ; 
    thread t2 {&Threading::consume,this} ; 
    t1.join() ; 
    t2.join() ; 
} 

void Threading::produce() { 
    random_device d;
    mt19937 mt(d());
    uniform_int_distribution<> distr(0., t_max);

    while (i<t_max) {
        this_thread::sleep_for(chrono::milliseconds(10));
        lock_guard<mutex> l(m); // lock_guard can only be locked and unlocked once. 
        int n = distr(mt);
        cout << "push: " << n << endl ; 
        q.push(n);
        i += n ;
        condition = true ; 
        cv.notify_all() ; // notify all condition_variables "cv" that are waiting that "condition = true"
    }
    cout << "all numbers pushed to queue" << endl ; ////// 
    while (!q.empty()) {
        condition = true ; 
        cv.notify_all() ;
    }
}


void Threading::consume() { 
    this_thread::sleep_for(chrono::milliseconds(100));
    unique_lock<mutex> l(m); // can be locked and unlocked multiple times

    while (i<t_max or !q.empty()) { // while the queue is not empty pop all numbers out of it
        cout << "pop:  " << q.front() << endl ; 
        int n = q.front() ;
        this_thread::sleep_for(chrono::milliseconds(n));
        q.pop() ;
        condition = false ;
        if (i<t_max or !q.empty()) {
            cv.wait(l, [&]{return condition ;}) ; // because condition is false "cv" is waiting until all numbers are pushed to queue
        }
    }
    cout << "all numbers popped from queue" << endl ; 
    // when the function is closed the mutex is released
}

#endif

namespace py = pybind11;
PYBIND11_MODULE(example, m) {
    py::class_<Threading>(m, "Threading")
        .def(py::init<double>(), "Constructor", py::arg("time") = 1.0) 
        .def("run", &Threading::run, "Run the code") ; //, "A class that runs two threads") 
}

