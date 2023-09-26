#include "String.hh"


String& String::operator=(const String& a) {
    if (this != &a) { insert(a.s) ;} // checks if not the string that is replaced by the other string is the same, otherwise this pointer will be deleted in the insert function. Then changes the old _s to the new _s inside "string"
    return *this ; // returns this string 
}


String& String::operator+=(const String& a) {
    int len_tot = len + a.len ; // set total length of new string
    char* str_tot = new char[len_tot] ; // creates pointer to total string (array of char)
    strcpy(str_tot,s) ; // set the first part of the new string to the value of the "s" variable that exists in the class String
    strcpy(str_tot+len,a.s) ; // set the second part of the new string (reached by setting pointer to: "str_tot + len") to the add string part 
    // also try the next line with {}:
    if (s) {delete[] s ;} // the old memory of s has to be delete. NOT: "if (s) {delete[] s}" otherwise the deleting takes place inside brackets
    s = str_tot ; // therefore the pointer s needs to be pointed to the new string, so that everything is still right in the class String
    len = len_tot ; // same for the length 
    return *this ;
} 


void String::insert(const char* str) { // private helper function
     len = strlen(str) ; // sets length of the input string
     if (s) delete[] s ; // deletes _s if it exists
     s = new char[len+1] ; // creates new _s of the right lenght
    strcpy(s,str) ; // copies old data to the new char pointer
  }

