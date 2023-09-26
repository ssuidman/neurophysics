#ifndef MANAGER 
#define MANAGER 

#include <string>
#include <iostream>
#include "Employee.hh"
using namespace std ;

class Manager : public Employee {
    public:
        Manager(const char* name, double salary,  set<Employee*> subs) : Employee(name,salary), subordinates(subs) {"Manager Constructor" ;}

        void addSubordinate(Employee& empl) {subordinates.insert(&empl) ;}
        const set<Employee*>& listOfSubordinates() const {return subordinates ;}
        void businessCard() const ;

    private:
        set<Employee*> subordinates ; 

} ;


void Manager::businessCard() const {
    Employee::businessCard() ;
    set<Employee*>::iterator iter ;
    cout << "   Managed employees: " ; 
    for (iter=subordinates.begin(); iter!=subordinates.end(); ++iter) {
        cout << (*iter)->name() << " " ; 
    }
    cout << endl ;
}

#endif 