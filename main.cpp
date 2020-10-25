#include <iostream>
#include "Matrix.h"
#include "matrix_utils.h"
#include "Function.h"
#include "list_4.h"
#include "list_5.h"
#include "interface.h"
#include "list_6.h"
#include "list_7.h"
#define ui unsigned int 

int main()
{
    
    int command;

    std::cout << "----- Welcome ------" << std::endl;
    
    list_commands();

    std::cout << "Enter a command: ";
    std::cin >> command;

    while(command){

        if(command == 1)
        {
            list_commands();
        }
        else if(command == 2)
        {
            decompose();
        }
        else if(command == 3)
        {
            determinant();
        }
        else if(command == 4)
        {
            solve_linear_system();
        }
        else if(command == 5)
        {
            eigen();

        }
        else if(command == 6)
        {
            mmq_call();
        }
        else if(command == 7)
        {
            test_positive_definite();
        }
        else if(command == 8)
        {
            root();
        }
        else if(command == 9)
        {
            solve_non_linear_equations();
        }
        else if(command == 10)
        {
            mmq_nl();
        }
        else if(command == 11)
        {
            integration();
        }
        else if(command == 12)
        {
            differentiate();
        }
        else if(command == 13)
        {
            ode();
        }
        list_commands();

        std::cout << "Enter a command: ";
        std::cin >> command;

    }

    return 0;
}