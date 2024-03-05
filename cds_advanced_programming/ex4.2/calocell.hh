#ifndef CALOCELL
#define CALOCELL

class Calocell {
    public:
        Calocell(double energy0, int ID0) {std::cout << "Calocell constructor " << this << std::endl ; energy = energy0; ID = ID0 ;}
        void set_energy(double new_energy) {energy = new_energy ;} 
        double get_energy() const {return energy ;} 
        int get_ID() const {return ID ;}
        void set_ID(int new_ID) {ID = new_ID ;}
    private:
        double energy = 1;
        int ID ;
} ;

#endif

