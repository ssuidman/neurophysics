#ifndef STRING_HH
#define STRING_HH
#include <cstring>

class String {
public:

  String(const char* str="") : s(0) { insert(str) ;}
  String(const String& a)  : s(0) { insert(a.s) ;}
  ~String() { delete[] s ;}
  String& operator=(const String& a) ; 
  String& operator+=(const String& a) ;
  operator const char*() const {return s ;} // Converts "const char*" to "String" and makes sure strlen() can be used. If this line were not there you get "error: cannot convert ‘String’ to ‘const char*’ "

  int length() const { return len ;}
  const char* data() const { return s ;}

private:

  char* s ;
  int len ;

  void insert(const char* str) ;

} ;

#endif 
