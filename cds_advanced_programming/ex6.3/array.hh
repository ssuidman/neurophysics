#ifndef ARRAY_HH
#define ARRAY_HH

template<class T>
class Array {
public:

  Array(int size, T default_val) : _size(size), _default(default_val) {_arr = new T[_size] ;} 
  Array(const Array& other) : _size(other._size), _default(other._default) {
    _arr = new T[other._size] ;
    for (int i=0 ; i<_size ; i++) { // Copy elements
      _arr[i] = other._arr[i] ;
    }
  }
  ~Array() {delete[] _arr ;}

  Array& operator=(const Array& other) 
  {
    if (&other==this) return *this ;
    if (_size != other._size) {
      resize(other._size) ;
    }
    for (int i=0 ; i<_size ; i++) {
      _arr[i] = other._arr[i] ;
    }
    return *this ;
  }

  T& operator[](int index) {
    if (index>=_size) {resize(index+1, _default) ;} // resize to index+1 if the index is bigger then the size. The "+1" comes from the fact that array[0] is also possible
    return _arr[index] ;
  }

  const T& operator[](int index) const {
    if (index>=_size) {resize(index+1, _default) ;}
    return _arr[index] ;
  }

  int size() const { return _size ; }

  void resize(int newSize, T default_val) {
    T* newArr = new T[newSize] ; // Allocate new array

    for (int i=0 ; i<_size ; i++) { // Copy elements
      newArr[i] = _arr[i] ;
    }
    for (int i = _size; i<newSize; i++) {
      newArr[i] = default_val ; 
    }

    delete[] _arr ; // Delete old array and install new one
    _size = newSize ;
    _arr = newArr ;
  }

private:
  int _size ;
  T* _arr ;
  T _default ;
} ;

#endif
