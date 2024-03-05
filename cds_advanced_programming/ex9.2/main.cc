#include <iostream>
#include <thread>
#include <queue>
#include <chrono>
#include <random>
#include <mutex>
#include <condition_variable>
using namespace std;

bool condition ; 
mutex m ; // the mutex can be owned by a thread. If the mutex is owned than all other threads are locked (so they are not continuing the code). The mutex can be passed over to different functions and there can be multiple mutex's
condition_variable cv ; // the condition variable is there to set a condition and is used often in combination with a lock that uses a mutex. 
int t_max = 3000;
int i = 0 ;

template<class T>
void produce(queue<T>& q) { 
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
        this_thread::sleep_for(chrono::milliseconds(n));
        condition = true ; 
        cv.notify_all() ; // notify all condition_variables "cv" that are waiting that "condition = true"
    }
    cout << "all numbers pushed to queue" << endl ; 
}

template<class T>
void consume(queue<T>& q) { 
    this_thread::sleep_for(chrono::milliseconds(100));
    unique_lock<mutex> l(m); // can be locked and unlocked multiple times

    while (i<t_max or !q.empty()) {
        cout << "pop:  " << q.front() << endl ; 
        int n = q.front() ;
        this_thread::sleep_for(chrono::milliseconds(n));
        q.pop() ;
        condition = false ;
        if (i<t_max or !q.empty()) {
            cv.wait(l, []{return condition ;}) ; // because condition is false "cv" is waiting until all numbers are pushed to queue
        }
    }
    cout << "all numbers popped from queue" << endl ; 
    // when the function is closed the mutex is released
}

int main() { 

    queue<double> Q ; 

    thread t1 {produce<double>,ref(Q)} ; 
    thread t2 {consume<double>,ref(Q)} ; 
    t1.join() ; 
    t2.join() ; 

    return 0 ; 
} 

// General:     Run the program with "g++ -o output main.cc -pthread", otherwise errors. 
// 9.2a)        You need to use std::ref() to get a reference passed in a thread. Passing 
//              by value is also easy. 

