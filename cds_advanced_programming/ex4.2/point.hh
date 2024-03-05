#ifndef POINT
#define POINT

class Point {

    public:
        Point() : x(0), y(0), z(0) {std::cout << "Point constructor " << this << std::endl ;}
        Point(double x0, double y0, double z0): x(x0), y(y0), z(z0) {std::cout << "Point constructor " << this << std::endl ;}
        double get_x() const {return x;}
        double get_y() const {return y;}
        double get_z() const {return z;}
        void set_x(double x_new) {x = x_new;}
        void set_y(double y_new) {y = y_new;}
        void set_z(double z_new) {z = z_new;}
        void set_xyz(double x_new, double y_new, double z_new) {x = x_new; y = y_new; z = z_new;}
    private:
        double x,y,z ;
} ;

#endif