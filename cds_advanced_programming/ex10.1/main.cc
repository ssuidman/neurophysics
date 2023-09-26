#include "functions.hh"

// All information regarding the functions is in functions.hh, but you need this main file to get it to bind to python. 
// The binding file can be created using: 
// g++ -O3 -Wall -shared -fPIC -std=c++11 -u dynamic_lookup `python3 -m pybind11  --includes` *.cc -o example`python3-config --extension-suffix`
// where "example" is the name of the module you are creating (import example). 
// Then you can start "python3" and run "import example", etc or just run the main.py file using "python3 main.py"
