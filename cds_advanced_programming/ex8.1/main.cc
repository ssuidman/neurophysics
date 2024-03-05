#include <iostream> 
#include <set> 
#include "Employee.hh" 
#include "Manager.hh" 
using namespace std ;

int main() {
    Employee empl("Sam 1",1000) ; // create employee 1
    empl.businessCard() ;
    cout << "   Salary: " << empl.salary() << endl << endl ; 

    Manager manager("Big Sam",2000,{&empl}) ; // Can also create empty set "{}" and then later use "addSubordate(empl)"
    
    Employee empl2("Sam 2", 1500) ; // create employee 2
    manager.addSubordinate(empl2) ; // give manager another subordinates empl2 

    set<Employee*> subs2 = manager.listOfSubordinates() ; // show the "listOfSubordinates" function
    (*subs2.find(&empl2))->businessCard() ; // The set subs2 contains references to employees. The "find" option makes sure you find the pointer of the set to a certain empl. To get the reference to an empl you need to use the "*" operator. You have now the reference(!) to an employee, so you need to use "->" to use a function inside the class "Employee". 
    cout << "   Salary: " << (*subs2.find(&empl2))->salary() << endl << endl ;
    

    manager.businessCard() ; // look at the manager again
    cout << "   Salary: " << manager.salary() << endl << endl ;

    Employee Wouter("Wouter",3000) ;
    Employee Ivo("Ivo",2000) ; 
    Manager Stan("Stan",2000,{}) ; // create managers without subordinates
    Manager Jo("Jo",1000,{}) ; 
    Manager Frank("Frank",1000,{}) ; 

    Stan.addSubordinate(Wouter) ; // add subordinates to the managers
    Stan.addSubordinate(Ivo) ;
    Frank.addSubordinate(Stan) ; // Even though Stan is a manager he is also an employee (inheritage). That is why he can easily be put in the list of subordinates of the (big) manager Frank 
    Frank.addSubordinate(Jo) ;
    
    Wouter.businessCard() ; 
    Ivo.businessCard() ; 
    Stan.businessCard() ; 
    Jo.businessCard() ; 
    Frank.businessCard() ; 


    return 0 ;

}



// 8.1c)    Better to use set in this over list/array, because there is no need for numbering 
//          your values and a set uses a hash table which has a faster lookup mechanism.
// 
// 8.1f)    It looks like the employee features of the manager behave exactly the same as for 
//          the other employees. 

