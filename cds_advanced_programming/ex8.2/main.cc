#include <iostream> 
#include <set> 
#include "Employee.hh" 
#include "Manager.hh" 
using namespace std ;

void printAllCards(set<Employee*>::const_iterator begin, set<Employee*>::const_iterator end) ;

int main() {
    Employee Wouter("Wouter",3000) ;
    Employee Ivo("Ivo",2000) ; 
    Manager Stan("Stan",2000,{}) ; // create managers without subordinates
    Manager Jo("Jo",1000,{}) ; 
    Manager Frank("Frank",1000,{}) ; 

    Stan.addSubordinate(Wouter) ; // add subordinates to the managers
    Stan.addSubordinate(Ivo) ;
    Frank.addSubordinate(Stan) ; // Even though Stan is a manager he is also an employee (inheritage). That is why he can easily be put in the list of subordinates of the (big) manager Frank 
    Frank.addSubordinate(Jo) ;
    set<Employee*> employees = {&Wouter, &Ivo, &Stan, &Jo, &Frank} ;
    printAllCards(employees.begin(),employees.end()) ;
    return 0 ;
}


void printAllCards(set<Employee*>::const_iterator begin, set<Employee*>::const_iterator end) {
    set<Employee*>::iterator iter ; 
    for(iter=begin; iter!=end; ++iter) {
        (*iter)->businessCard() ;
    }
}

// 8.2d)    The output is now different, because the businesscards are printed the 
//          same as if there are no managers (the subordinates are not printed). 
// 
// 8.2g)    After adding "virtual" the code behaves differently and prints the 
//          manager cards as manager cards even though we did not specify they were
//          managers. 