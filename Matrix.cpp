#include "Matrix.h"
#include "matrix_utils.h"
#include "interface.h"
#include <vector>

using vectors = std::vector<double>;
using matrix = std::vector<vectors>;

unsigned int Matrix::get_n()
{
    return this->n;
}

unsigned int Matrix::get_m()
{
    return this->m;
}

std::string Matrix::get_name()
{
    return this->name;
}

double Matrix::get_determinant()
{
    return this->determinant;
}

void Matrix::set_n(unsigned int n)
{
    this->n = n;    
}

void Matrix::set_m(unsigned int m)
{
    this->m = m;
}

void Matrix::set_name(std::string name)
{
    this->name = name;
}

void Matrix::set_determinant(double determinant)
{
    this->determinant = determinant;
}


Matrix::Matrix(unsigned int n, unsigned int m, std::string name)
{
    this->set_n(n);
    this->set_m(m);
    this->set_name(name);
    this->entries = matrix(this->get_n(), vectors(this->get_m()));
}

void Matrix::read()
{
    divisor();

    std::cout << "Type the entries of the Matrix from top-left to bottom right: " << std::endl;

    for(int i = 0; i < this->get_n(); i++)
    {
        for(int j = 0; j<this->get_m(); j++)
        {
            std::cin >> this->entries[i][j];
        }
    }

    divisor();
}

void Matrix::write()
{
    divisor();

    std::cout << '[' << this->get_name() << ']' << std::endl;
    
    for(int i = 0; i < this->get_n(); i++)
    {
        for(int j = 0; j < this->get_m(); j++)
        {
            double aij = this->entries[i][j];
            printf("%.4lf ",aij);
        }
        std::cout << std::endl;
    }
    
    divisor();
}

void Matrix::id()
{
    for(int i = 0; i < this->get_n(); i++)
    {
        for(int j = 0; j < this->get_m(); j++)
        {
            if(i == j){
                this->entries[i][j] = 1;
            }
            else
            {
                this->entries[i][j] = 0;
            }
        }
    }
}


Matrix::~Matrix()
{
    
}