#include <vector>
#include <iostream>

using vectors = std::vector<double>;
using matrix = std::vector<vectors>;

class Matrix
{
    public:
        Matrix                      (unsigned int n,unsigned int m,std::string name);
        ~Matrix                     ();
        
        matrix entries;

        unsigned int get_n          ();
        unsigned int get_m          ();
        std::string get_name        ();
        double get_determinant         ();
        
        void set_n                  (unsigned int n);
        void set_m                  (unsigned int m);
        void set_name               (std::string name);
        void set_determinant        (double determinant);

        void write                  ();
        void read                   (); 
        void id                     ();

    private:
        unsigned int n;
        unsigned int m;
        std::string name;
        double determinant;
        

};