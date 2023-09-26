#include <iostream> 
#include "Telephone.hh" 

int main() {
    Telephone telephone ; 
    Telephone telephone2(telephone) ; // copy constructor used
    return 0 ;
}

// 4.1a)    See the file "tree_diagram.txt"
// 
// 4.1b)    you need to put this line here, because otherwise in "Telephone.hh" "Cable.hh" is 
//          included (and thus the class is defined here) and again when calling "Handset.hh" 
//          there is "Cable.hh" included and again defined the class. That is why you need to 
//          set these lines here so that the second definition of the class is skipped
// 
// 4.1c)    When calling the class "Telephone" first classes at the bottom of the tree are called. 
//          Then the level above this all the way up to the Telephone class. This is because the 
//          classes that use another class first need to get instances of the other classes to 
//          create their own instances and that can only be done if the other classes are already 
//          are activated. This can be seen in the way constructors are called when running this 
//          line: Cable-->Chassis-->Shell-->Housing-->Button[12]-->etc. The destructor has an 
//          exactly reversed order. This makes sense, because if at the bottom of a tree something 
//          is destructed this means that the upper classes that are not yet destructed cannot use 
//          this class anymore. This is unwanted of course. 
// 
// 4.1d)    Copy constructor used above in main.cc
// 
// 4.1e)    The instances that are created (buttton for example) to create the copied telephone 
//          object "telephone2" are not copied ones from the previous instances. This means that a 
//          new button is created with its own properties. You want to copy a full object with the 
//          same properties that is why you need to use the copy constructors for all the objects. 
// 
// 4.2f)    It works after using " ...: cable(telephone.cable), ... "
// 
// 4.2g)    Dialer class is changed with new operator and delete operator to not copy the pointer. 
//          If you would copy the pointer, changing memory in one instance would also change things 
//          in the other(s).
// 
// 4.2h)    Constructor/copy constructor/destructor all changed. 
