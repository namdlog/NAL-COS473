#include <iostream>
#include "Matrix.h"
#include <string>
#include <math.h>
#include "matrix_utils.h"
#include "Function.h"
#include "list_4.h"
#include "list_5.h"
#include "interface.h"
#include "list_6.h"
#include "list_7.h"
#include <iomanip>      

Matrix enter(std::string name)
{
    
    int n;
    std::cout << "Enter the first dimention of the matrix (number of rows): ";
    std::cin >> n;
    int m;
    std::cout << "Enter the second dimention of the matrix (number of columns): ";
    std::cin >> m;
    
    Matrix A(n,m,name);

    A.read();

    A.write();
    
    return A;
}

Function new_function(std::string name)
{
    
    unsigned int degree;
    std::cout << "Enter the degree of the function: ";
    std::cin >> degree;
    
    Function f(degree,name);

    f.read();
    
    return f;
}

void divisor()
{
    std::cout << "--------------------" << std::endl;
}

void list_commands()
{
    divisor();

    std::cout << " 0 - Exit\n 1 - List commands\n 2 - Decompose\n 3 - Calculate the determinant\n 4 - Solve linear system\n 5 - Eigenvalues and Eigenvectors\n 6 - MMQ\n 7 - Test Definite Positive\n 8 - Calculate the Root of a Function\n 9 - Solve non linear system\n 10 - MMQ NL\n 11 - Integration\n 12 - Differentiate\n 13 - ODE" << std:: endl;

    divisor();
}

void solve_linear_system()
{
    int option;

    Matrix A = enter("A");

    Matrix b = enter("b");

    int n = b.get_n();
    int m = A.get_m();

    Matrix x(n,1,"x");

    divisor();

    std::cout << " 1 - LU\n 2 - Cholesky\n 3 - Ax = b\n 4 - Iterative: Jacobi\n 5 - Iterative: Gauss-Seidel" << std:: endl;
    
    divisor();
    
    std::cout << "Choose a method: ";
    std::cin >> option;

    if(option == 1)
    {

        Matrix L(n,m,"L");
        Matrix U(n,m,"U");
        divisor();

        if(decomposition_lu(A,L,U))
        {
            L.write();
            U.write();
            solve_linear_system_trick(L,U,x,b);
            U.write();
            x.write();
        }
        else
        {   
            std::cout << "Sorry, try again" << std::endl;   
        }

        divisor();

    }
    else if(option == 2)
    {   
        Matrix L(n,m,"L");
        Matrix Lt(n,m,"Lt");
        
        divisor();

        if(decomposition_cholesky(A,L,Lt))
        {
            solve_linear_system_trick(L,Lt,x,b);
            x.write();
        }
        else
        {   
            std::cout << "Sorry, try again" << std::endl;   
        }

        divisor();

    }else if(option == 3)
    {    

        Matrix M(n,m,"M");
        M.id();

        Matrix U(n,m,"U");

        gauss_elimination(A,M,U);
        

        multiply(M,b,b);
    
        
        if(solve_linear_system(U,x,b))
        {
            x.write();
        }else
        {
            std::cout << "Sorry, try again" << std::endl;
        }

    }else if(option == 4)
    {
        double tol;
        
        std::cout << "Enter the tolerante: " << std::endl;
        std::cin >> tol;

        jacobi(A,x,b,tol);
 
        x.write();

    }else if(option == 5)
    {
        double tol;
        
        std::cout << "Enter the tolerante: " << std::endl;
        std::cin >> tol;

        gauss_seidel(A,x,b,tol);
 
        x.write();
 
    }

    divisor();
}

void decompose()
{
    int option;

    Matrix A = enter("A");

    int n = A.get_n();
    int m = A.get_m();

    divisor();

    std::cout << " 1 - LU\n 2 - Cholesky" << std:: endl;
    
    divisor();
    
    std::cout << "Choose a method: ";
    std::cin >> option;

    if(option == 1)
    {

        Matrix L(n,m,"L");
        Matrix U(n,m,"U");
        
        divisor();

        if(decomposition_lu(A,L,U))
        {
            L.write();
            U.write();
        }
        else
        {   
            std::cout << "Sorry, try again" << std::endl;   
        }

        divisor();
    }
    else if(option == 2)
    {   
        Matrix L(n,m,"L");
        Matrix Lt(n,m,"Lt");
        if(decomposition_cholesky(A,L,Lt))
        {
            L.write();
            Lt.write();
        }
        else
        {   
            std::cout << "Sorry, try again" << std::endl;   
        }
    }

    divisor();
}

void determinant()
{
    Matrix A = enter("A");

    int n = A.get_n();
    int m = A.get_m();

    divisor();

    A.set_determinant(calculate_determinant(A));
    printf("The determinant of the Matrix is: %.8lf", A.get_determinant());
    std::cout << std::endl;
    
    
    divisor();
    
}

void eigen()
{
    int option;
    Matrix A = enter("A");

    int n = A.get_n();
    int m = A.get_m();

    divisor();

    std::cout << " 1 - Power Method\n 2 - Jacobi" << std:: endl;
    
    std::cout << "Choose a method: ";
    std::cin >> option;

    if(option == 1)
    {
        Matrix xz(n,1,"Eigenvector");
        
        double max_eigen_value = 1.0;
        double tol;
        
        std::cout << "Enter the tolerance: " << std::endl;
        std::cin >> tol;

        if(eigen_power_method(A, xz, max_eigen_value, tol))
        {
            xz.write();
            std::cout << "The greatest eigen value of the matrix is: " << max_eigen_value << std::endl;
        }else
        {
            std:: cout << "Sorry, try again" << std::endl;
        }
    }else if(option == 2)
    {
        Matrix Ak(n,n,"Eigenvalues");
        Matrix Xk(n,n,"Eigenvectors");
        Xk.id();
        double tol;
        std::cout << "Enter the tolerance: " << std::endl;
        std::cin >> tol;

        if(eigen_jacobi(A, Ak, Xk,tol))
        {
            Ak.write();
            Xk.write();
        }
    }
    divisor();

}

void mmq_call()
{
    divisor();  
    
    std::vector<std::pair<double,double>> samples;
    std::vector<std::vector<double>> coeficients;
    
    int n,m;
    double x_p,y_p,phi;
    
    std::cout << "Enter the number of sample points: ";
    std::cin >> n;
    std::cout << "Enter the x and y coordenates: " << std::endl;
    for(int i=0;i < n; i++){
        std::cout << "X: ";
        std::cin >> x_p;
        std::cout << "Y: ";
        std::cin>> y_p;
        samples.push_back({x_p,y_p});
    }
    
    bool linear = false;
    std::string answer;
    std::cout << "Do you want linear? Yes or No: ";
    std::cin >> answer;
    if(answer == "Yes")
    {
        linear = true;
    }

    if(!linear)
    {
        std::cout << "Enter the number of coeficients: ";
        std::cin >> m;
        Matrix x(m,1,"Function");
        std::cout << "Enter the x and y coordenates: " << std::endl;
        for(int i = 0; i < n; i++)
        {
            std::vector<double> coef;
            for(int j = 0; j < m; j++)
            {
                std::cout << "phi" << i << "(x" << j << "): ";
                std::cin >> phi;
                coef.push_back(phi);
            }
            coeficients.push_back(coef);
        }
        mmq_calc(x,samples,coeficients);

        x.write();
    }else
    {
        Matrix x(2,1,"Function");
        for(int i = 0; i < n; i++)
        {
            std::vector<double> coef;
            coef.push_back(1);
            coef.push_back(samples[i].first);
            coeficients.push_back(coef);
        }
        mmq_calc(x,samples,coeficients);

        x.write();
    }
   
    
    divisor();
}

void test_positive_definite()
{
    divisor();
    
    Matrix A = enter("A");

    if(positive_definite(A))
    {
        std::cout << "The matrix is positive definite" << std::endl;
    }else
    {  
        std::cout << "The matrix is not positive definite" << std::endl;
    }

    divisor();

}

double function_slide1(double x){
    return pow(x,2)-4*cos(x);
}

double function_l4q1(double x){
    double g = 9.806;
    double k = 0.00341;
    return log(cosh(x*sqrt(g*k)))-50.0;
}

double function_l4q2(double x){
    return 4*cos(x)-exp(2*x);
}


double function_exec(double x){
    return function_l4q1(x);
}

void root()
{
    int option;

    divisor();
    std::cout << "Wich method do you want to use: " << std::endl; 

    std::cout << " 1 - Bissection\n 2 - Newton\n 3 - Secant\n 4 - Inverse Interpolation" << std::endl;
    std::cin >> option;

    if(option == 1)
    {
        double a,b,tol;
    
        std::cout << "Enter the limit a: ";
        std::cin >> a;
        std::cout << "Enter the limit b: ";
        std::cin >> b;
        std::cout << "Enter the tolerance: ";
        std::cin >> tol;

        bisectionMethod(function_exec,a,b,tol);
    }
    else if(option == 2)
    {
        double xZero,tol;
        std::cout << "Enter with xZero :";
        std::cin >> xZero;
        std::cout << "Enter with the tolerance: ";
        std::cin >> tol;
        
        newtonMethod(function_exec,xZero,tol,1000);
    }
    else if(option == 3)
    {
        double xZero,tol;
        std::cout << "Enter with xZero :";
        std::cin >> xZero;
        std::cout << "Enter with the tolerance: ";
        std::cin >> tol;
        
        secantMethod(function_exec,xZero,tol,1000);

    }else if(option == 4)
    {
        int nSamples,inter;
        double tol;

        std::cout << "Enter the number of points: ";
        std::cin >> nSamples;
        
        std::vector<double> x(nSamples);
        
        for(int i=0;i<nSamples;i++){
            std::cout << "Enter x[" << i << "]: ";
            std::cin >> x[i];
        }
        
        std::cout << "Enter the tolerance: ";
        std::cin >> tol;
        
        inter = 10000;
        rootInverseInterpolation(function_exec,x,tol,inter);
    }
    /*Function f = new_function("f");
    
    f.write();
    
    double x;

    std::cout << "X pra calcular :" << std::endl;
    std::cin >> x;
    std::cout << f.evaluation(x) << std::endl;*/
    
    divisor();
}

double nl_1(int offset, vectors variables)
{
    switch(offset)
    {
        case 0:
            return 16.0*pow(variables[0],4)+16.0*pow(variables[1],4)+pow(variables[2],4)-16.0;
        case 1:
            return pow(variables[0],2)+pow(variables[1],2)+pow(variables[2],2)-3.0;
        case 2:
            return pow(variables[0],3)-variables[1]+variables[2]-1.0;
    }
}

double nl_2(int offset, vectors variables)
{
    switch(offset)
    {
        case 0:
            return variables[0]+2*variables[1]-2.0;
        case 1:
            return pow(variables[0],2)+4*pow(variables[1],2)-4.0;
    }
}

double nl_3(int offset,vectors variables)
{   
    double theta_1 = 0.0;
    double theta_2 = 3.0;

    switch(offset)
    {
        case 0:
            return 2*pow(variables[1],2)+pow(variables[0],2)+6*pow(variables[2],2)-1.0;
        case 1:
            return 8*pow(variables[1],3)+6*variables[1]*pow(variables[0],2)+36*variables[1]*variables[0]*variables[2]+108*variables[1]*pow(variables[2],2)-theta_1;
        case 2:
            return 60*pow(variables[1],4)+60*pow(variables[1],2)*pow(variables[0],2)+576*pow(variables[1],2)*variables[0]*variables[2]+2232*pow(variables[1],2)*pow(variables[2],2)+252*pow(variables[2],2)*pow(variables[0],2)+1296*pow(variables[2],2)*variables[0]+3348*pow(variables[2],4)*24*pow(variables[0],3)*variables[2]+3*variables[0]-theta_2;
    }
}

double nl_function(int offset,  vectors variables)
{
    return nl_3(offset,variables);
}

void solve_non_linear_equations()
{
    divisor();

    int n_equations = 3;
    int m_variables = 3;
    int iterationNumber;
    double tol;
    
    Matrix x_zero(m_variables,1,"x_zero");
    Matrix x_answer(m_variables,1,"x_answer");
    
    for(int i = 0; i < m_variables; i++)
    {
    
        std::cout << "Enter with the value of x[" << i << "]: ";
        std::cin >> x_zero.entries[i][0];
    
    
    }
    x_zero.write();
    std::cout << "Enter the tolerance: ";
    std::cin >> tol;

    iterationNumber = 10000;

    std::cout << "Wich method do you want? 1 - Newton 2 - Broyden" << std::endl;
    int method;
    std::cin >> method;
    
    if(method == 1)
    {
        x_answer = newtonsMethodNlq(n_equations,nl_function,x_zero,tol,iterationNumber);
        x_answer.write();
    }
    else
    {
        x_answer = broydenMethodNlq(n_equations,nl_function,x_zero,tol,iterationNumber);
        x_answer.write();
    }
    
    divisor();
}

double function_l4q5(int offset, vectors variables){
    
    if(offset == 0){
        std::pair<double,double> point = {1.0, 1.0};
        double x = point.first;
        double y = point.second;
        return variables[0] + variables[1]*pow(x,variables[2]) - y;
    }
    if(offset == 1){
        std::pair<double,double> point = {2.0, 2.0};
        double x = point.first;
        double y = point.second;
        return variables[0] + variables[1]*pow(x,variables[2]) - y;
    }
    if(offset == 2){
        std::pair<double,double> point = {3.0, 9.0};
        double x = point.first;
        double y = point.second;
        return variables[0] + variables[1]*pow(x,variables[2]) - y;
    }

}

void mmq_nl()
{
    //Matrix mmqNl(int n_equations, double function(int offset, vectors variables), Matrix x_zero,double tolerance,int iterationsNumber)
    
    Matrix x(3,1,"x");
    Matrix xzero(3,1,"xZero");
    
    std::cout << "Enter with the x zero values: " << std::endl;
    std::cin >> xzero.entries[0][0];
    std::cin >> xzero.entries[1][0];
    std::cin >> xzero.entries[2][0];

    x = mmqNl(3,function_l4q5,xzero,0.0001,1000);
    x.write();
}

double PI = 3.1415926535897932384626433832795028841971693993751058209;

double function_slide2(double x)
{
    return exp(-1.0*pow(x,2));
}

double function_l5q2(double x)
{
    return exp(pow(x,2)/-2.0)*(1/sqrt(2*PI));
}

double function_l5q3(double x)
{
    double omegn = 1.0;
    double eps   = 0.05;
    return 2.0*pow(1.0/sqrt(pow(1-pow(x/omegn,2),2)+pow(2*eps*x/omegn,2)),2);
}

double function_l5q3_2(double x)
{
    double omegn = 1.0;
    double eps   = 0.05;
    return pow(x,2)*2.0*pow(1.0/sqrt(pow(1-pow(x/omegn,2),2)+pow(2*eps*x/omegn,2)),2);
}

double function_l5q4(double x)
{
    double omegn = 1.0;
    double eps   = 0.05;
    double Hs    = 3.0;
    double Tz    = 5.0;

    double sn = ((4*pow(PI,3)*pow(Hs,2))/(pow(x,5)*pow(Tz,4)))*exp((-16.0*pow(PI,3))/(pow(x,4)*pow(Tz,4)));
    double ss = pow(1.0/sqrt(pow(1-pow(x/omegn,2),2)+pow(2*eps*x/omegn,2)),2);
    return ss*sn;
}

double function_l5q4_2(double x)
{
    return x*x*function_l5q4(x);
}

double function_l5q5(double x)
{
    return 2+2*x-(x*x)+3*pow(x,3);
}

double function_l5q6(double x)
{
    return 1.0/(1.0+pow(x,2));
}

double function_l52_q1_i1(double x)
{
    return pow(x,3)+1.0/exp(x);
}

double function_l52_q1_i2(double x)
{
    return pow(x,1/3)+log(x);
}

double function_l52_q1_i3(double x)
{  
    return 1-exp(-1.0*(pow(x/5,2)));
}

double function_derivate(double x)
{
    return function_l52_q1_i3(x);
}

std::vector<double> weight_gauss[11][2],weight_hermite[11][2],weight_laguerre[11][2];
    
void iniciate_gauss(std::vector<double> (&W)[11][2]){
    //Gauss
    W[2][0] = {1.0000000000000000,1.0000000000000000};
    W[2][1] = {-0.5773502691896257,0.5773502691896257}; 
    W[3][0] = {0.8888888888888888,0.5555555555555556,0.5555555555555556};
    W[3][1] = {0.0,-0.7745966692414834,0.7745966692414834};
    W[4][0] = {0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538};
    W[4][1] = {-0.3399810435848563,0.3399810435848563,-0.8611363115940526,0.8611363115940526};
    W[5][0] = {0.5688888888888889,0.4786286704993665,0.4786286704993665,0.2369268850561891,0.2369268850561891};
    W[5][1] = {0.0000000000000000,-0.5384693101056831,0.5384693101056831,-0.9061798459386640,0.9061798459386640};   
    W[6][0] = {0.3607615730481386,0.3607615730481386,0.4679139345726910,0.4679139345726910,0.1713244923791704,0.1713244923791704};
    W[6][1] = {0.6612093864662645,-0.6612093864662645,-0.2386191860831969,0.2386191860831969,-0.9324695142031521,0.9324695142031521};
    W[7][0] = {0.4179591836734694,0.3818300505051189,0.3818300505051189,0.2797053914892766,0.2797053914892766,0.1294849661688697,0.1294849661688697};
    W[7][1] = {0.0,0.4058451513773972,-0.4058451513773972,-0.7415311855993945,0.7415311855993945,-0.9491079123427585,0.9491079123427585};
    W[8][0] = {0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763};
    W[8][1] = {-0.1834346424956498,0.1834346424956498,-0.5255324099163290,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363};
    W[9][0] = {0.3302393550012598,0.1806481606948574,0.1806481606948574,0.0812743883615744,0.0812743883615744,0.3123470770400029,0.3123470770400029,0.2606106964029354,0.2606106964029354};
    W[9][1] = {0.0,-0.8360311073266358,0.8360311073266358,-0.9681602395076261,0.9681602395076261,-0.3242534234038089,0.3242534234038089,-0.6133714327005904,0.6133714327005904};
    W[10][0] = {0.2955242247147529,0.2955242247147529,0.2692667193099963,0.2692667193099963,0.2190863625159820,0.2190863625159820,0.1494513491505806,0.1494513491505806,0.0666713443086881,0.0666713443086881};
    W[10][1] = {-0.1488743389816312,0.1488743389816312,-0.4333953941292472,0.4333953941292472,-0.6794095682990244,0.6794095682990244,-0.8650633666889845,0.8650633666889845,-0.9739065285171717,0.9739065285171717};
    
}

void iniciate_gauss_hermite(std::vector<double> (&W)[11][2]){
   W[2][1] = {-0.7071067811865475244008,0.7071067811865475244008};
   W[2][0] = {0.8862269254527580136491,0.8862269254527580136491};
   W[3][1] = {-1.224744871391589049099,0.0,1.224744871391589049099};
   W[3][0] = {0.295408975150919337883,1.181635900603677351532,0.295408975150919337883};
   W[4][1] = {-1.650680123885784555883,-0.5246476232752903178841,0.5246476232752903178841,1.650680123885784555883};
   W[4][0] = {0.081312835447245177143,0.8049140900055128365061,0.8049140900055128365061,0.08131283544724517714303};
   W[5][1] = {2.020182870456085632929,-0.9585724646138185071128,0.0,0.9585724646138185071128,2.020182870456085632929};
   W[5][0] = {0.01995324205904591320774,0.3936193231522411598285,0.9453087204829418812257,0.393619323152241159828,0.01995324205904591320774};
   W[6][1] = {-2.350604973674492222834,-1.335849074013696949715,-0.4360774119276165086792,0.436077411927616508679,1.335849074013696949715,2.350604973674492222834};
   W[6][0] = {0.004530009905508845640857,0.1570673203228566439163,0.7246295952243925240919,0.724629595224392524092,0.1570673203228566439163,0.004530009905508845640857};
   W[7][1] = {-2.350604973674492222834,-1.335849074013696949715,-0.4360774119276165086792,0.436077411927616508679,1.335849074013696949715,2.350604973674492222834};
   W[7][0] = {0.004530009905508845640857,0.1570673203228566439163,0.7246295952243925240919,0.724629595224392524092,0.1570673203228566439163,0.004530009905508845640857},
   W[8][1] = {-2.930637420257244019224,-1.981656756695842925855,-1.157193712446780194721,-0.3811869902073221168547,0.3811869902073221168547,1.157193712446780194721,1.981656756695842925855,2.930637420257244019224};
   W[8][0] = {1.99604072211367619206E-4,0.0170779830074134754562,0.2078023258148918795433,0.6611470125582412910304,0.6611470125582412910304,0.2078023258148918795433,0.0170779830074134754562,1.996040722113676192061E-4},
   W[9][1] = {-3.19099320178152760723,-2.266580584531843111802,-1.468553289216667931667,-0.7235510187528375733226,0.0,0.7235510187528375733226,1.468553289216667931667,2.266580584531843111802,3.19099320178152760723};
   W[9][0] = {3.960697726326438190459E-5,0.00494362427553694721722,0.088474527394376573288,0.4326515590025557501998,0.7202352156060509571243,0.4326515590025557502,0.088474527394376573288,0.004943624275536947217225,3.96069772632643819046E-5},
   W[10][1] = {-3.436159118837737603327,-2.532731674232789796409,-1.756683649299881773451,-1.036610829789513654178,-0.3429013272237046087892,0.3429013272237046087892,1.036610829789513654178,1.756683649299881773451,2.532731674232789796409,3.436159118837737603327};
   W[10][0] = {7.64043285523262062916E-6,0.001343645746781,0.03387439445548,0.2401386110823146864165,0.61086263373257987836,0.6108263373532,0.24013861108,0.033874394455481,0.00134364574812326,7.64043285523262062916E-6};
}

void iniciate_laguerre(std::vector<double> (&W)[11][2]){
    W[2][1] = {0.5857864376269049511983,3.414213562373095048802};
    W[2][0] = {0.8535533905932737622004,0.1464466094067262377996};
    W[3][1] = {0.4157745567834790833115,2.29428036027904171982,6.289945082937479196866};
    W[3][0] = {0.71109300992917301545,0.2785177335692408488014,0.01038925650158613574897};
    W[4][1] = {0.3225476896193923118004,1.745761101158346575687,4.536620296921127983279,9.395070912301133129234};
    W[4][0] = {0.603154104341633601636,0.3574186924377996866415,0.03888790851500538427244,5.392947055613274501038E-4};
    W[5][1] = {0.2635603197181409102031,1.413403059106516792218,3.596425771040722081223,7.085810005858837556922,12.64080084427578265943};
    W[5][0] = {0.5217556105828086524759,0.3986668110831759274541,0.0759424496817075953877,0.00361175867992204845446,2.33699723857762278911E-5};
    W[6][1] = {0.2228466041792606894644,1.188932101672623030743,2.992736326059314077691,5.77514356910451050184,9.837467418382589917716,15.98287398060170178255};
    W[6][0] = {0.4589646739499635935683,0.417000830772120994113,0.113373382074044975739,0.01039919745314907489891,2.610172028149320594792E-4,8.9854790642962123883E-7};
    W[7][1] = {0.1930436765603624138383,1.026664895339191950345,2.567876744950746206908,4.900353084526484568102,8.182153444562860791082,12.73418029179781375801,19.39572786226254031171};
    W[7][0] = {0.4093189517012739021304,0.4218312778617197799293,0.1471263486575052783954,0.02063351446871693986571,0.001074010143280745522132,1.586546434856420126873E-5,3.17031547899558056227E-8};
    W[8][1] = {3.17031547899558056227E-8,0.903701776799379912186,2.25108662986613068931,4.266700170287658793649,7.045905402393465697279,10.75851601018099522406,15.74067864127800457803,22.8631317368892641057};
    W[8][0] = {0.369188589341637529921,0.418786780814342956077,0.1757949866371718056997,0.03334349226121565152213,0.00279453623522567252494,9.07650877335821310424E-5,8.48574671627253154487E-7,1.04800117487151038162E-9};
    W[9][1] = {0.152322227731808247428,0.8072200227422558477414,2.005135155619347122983,3.783473973331232991675,6.204956777876612606974,9.37298525168757620181,13.4662369110920935711,18.83359778899169661415,26.37407189092737679614};
    W[9][0] = {0.3361264217979625196735,0.4112139804239843873091,0.1992875253708855808606,0.0474605627656515992621,0.005599626610794583177004,3.05249767093210566305E-4,6.59212302607535239226E-6,4.1107693303495484429E-8,3.29087403035070757647E-11};
    W[10][1] = {0.1377934705404924308308,0.72945454950317049816,1.8083429017,3.40143369785,5.552496,8.330152746,11.843785837900,16.27925783137810209953,21.99658581198076195128,29.92069701227389155991};
    W[10][0] = {0.3084411157650201415475,0.401119929155273551516,0.2180682876118094215886,0.0620874560986777473929,0.00950151697518110055384,7.53008388587538775456E-4,2.82592334959956556742E-5,4.24931398496268637259E-7,1.839564823979630780922E-9,9.911827219609008558378E-13};
}

double function_l5q7_1(double x)
{
    return (1/sqrt(2*PI))*exp(-0.5*pow(x,2));
}

double function_l5q7_2(double x)
{
    return exp(-0.5*pow(x,2))*(pow(x,2)/sqrt(2*PI));
}

double function_integration(double x)
{
    return function_l5q7_1(x);
}

double laguerre_function(double x)
{
    return function_integration(x)/exp(-1.0*x);
}

double hermite_function(double x)
{
    double a = function_integration(x);
    double b = exp(-1*pow(x,2));
    //std::cout << "X: " << x << std::endl;
    //std::cout << "Resultado = " << a << " " << b << " " << a/b << std::endl;
    return a/b;
}

void integration()
{

    int nPoints;
    int method;
    std::string a_limit,b_limit;
    double hermite_integral = 0.0;
    double laguerre_integral = 0.0;
    double a;
    double b;
    double answer;
    
    std::cout << "Please enter the number of points: " << std::endl;
    std::cin >> nPoints;
    std::cout << "Choose your method: 1 - Polynomial 2 - Gauss" << std::endl;
    std::cin >> method;
    
    if(method == 2){
        iniciate_gauss(weight_gauss);
        iniciate_gauss_hermite(weight_hermite);
        iniciate_laguerre(weight_laguerre);
    }

    std::cout << "Enter with the limits: " << std::endl;
    std::cout << "a: ";
    std::cin >> a_limit;
    std::cout << "b: ";
    std::cin >> b_limit;

   
    for(int i=0;i<nPoints;i++){
        laguerre_integral += laguerre_function(weight_laguerre[nPoints][1][i])*weight_laguerre[nPoints][0][i];
    }
    for(int i=0;i<nPoints;i++){
        //std::cout << weight_hermite[nPoints][1][i] << " " << weight_hermite[nPoints][0][i] << std::endl;
        //std::cout << hermite_integral << std::endl;
        //std::cout << "Func Test: " << function_integration(weight_hermite[nPoints][1][i]) << std::endl;
        //std::cout << "Func Off: " << hermite_function(weight_hermite[nPoints][1][i]) << std::endl;
        
        //std::cout << "Hermite acumulado " << hermite_integral << std::endl;
        hermite_integral  += hermite_function(weight_hermite[nPoints][1][i])*weight_hermite[nPoints][0][i];
    }
   
    if(a_limit == "-inf" or b_limit == "inf"){
        if(a_limit == "-inf"){
            if(b_limit == "inf"){
                answer = hermite_integral;
            }else{
                b = stoi(b_limit);
                if(b>=0){
                    std::cout << "Entrei aqui" << std::endl;
                    std::cout << "Hermite: " << hermite_integral << std::endl;
                    std::cout << "Laguerre: " << laguerre_integral << std::endl;
                    answer = hermite_integral - laguerre_integral + gaussIntegration(function_integration,0,b,nPoints,weight_gauss);
                }else{             
                    answer = laguerre_integral;
                }
            }
        }
        else{
            a = stoi(a_limit);
            if(a<0){
                answer = laguerre_integral + gaussIntegration(function_integration,a,0,nPoints,weight_gauss);
            }else if(a>0){
                answer = laguerre_integral - gaussIntegration(function_integration,0,a,nPoints,weight_gauss);
            }else{
                answer = laguerre_integral;
            }
        }
    }else{
        a = stoi(a_limit);
        b = stoi(b_limit);
        if(method == 1)
        {
            answer = polynomialIntegration(function_integration,a,b,nPoints);
        }
        else
        {
            answer = gaussIntegration(function_integration,a,b,nPoints,weight_gauss);
        }
    }
    std::cout << "The integration is: " << answer << std::endl;
}

void differentiate()
{

    divisor();

    bool method = false;
    int type = 0;
    int p = 0;
    double x;
    double ans = 0.0;
    std::cout << "Do you want to use Richard? 0 - No / 1 - Yes" << std::endl;
    std::cin >> method;
    std::cout << "Wich method do you want? 0 - Back / 1 - Central / 2 - Foward" << std::endl;
    std::cin >> type;
    if(method)
    {
        std::cout << "What is the value of p: " << std::endl;
        std::cin >> p;
    }
    std::cout << "Enter the value of x: " << std::endl;
    std::cin >> x;
    double delta_x = 0.5;
    if(method)
    {
        if(type == 0)
        {
            ans = derivadaPassoTrasRichard(function_derivate,x,delta_x,p);
        }else if(type == 1)
        {
            ans = derivadaPassoCentralRichard(function_derivate,x,delta_x,p);
        }else if(type == 2)
        {
            ans = derivadaPassoFrenteRichard(function_derivate,x,delta_x,p);
        }
    }
    else
    {
        if(type == 0)
        {
            ans = derivadaPassoTras(function_derivate,x,delta_x);
        }
        else if(type == 1)
        {
            ans = derivadaPassoCentral(function_derivate,x,delta_x);
        }
        else if(type == 2)
        {
            ans = derivadaPassoFrente(function_derivate,x,delta_x);
        }
    }

    std::cout << ans << std::endl;


    divisor();
}

double l6_q1(double t,double y)
{
    return -2*t*(pow(y,2));
}
double diferential_function(double t,double y)
{
        return l6_q1(t,y);
}

double l7_q2(double t,double y,double yl)
{
    return 2*sin(0.5*t)+sin(t)+cos(1.5*t)-0.2*yl-y;
}

double l7_q3(double t,double y,double yl)
{
    return -9.81 - yl*pow(pow(yl,2),0.5);
}

double diferential_function_second(double t,double y, double yl)
{
        return l7_q2(t,y,yl);
}

void ode()
{
    divisor();
    
    int method;
    double deltaT;
    std::pair<std::vector<double>,std::vector<double>> ans;
    std::cout << "Wich method do you want to use? 0 - Euler / 1 - Runge-Kutta Second Order / 2 - Runge-Kutta - Forth Order / 3 - Taylor / 4 - Runge Kutta Nystron" << std::endl;
    std::cin >> method;
    std::cout << "Enter with delta: ";
    std::cin >> deltaT;
    double ini              = 0.0;
    double final            = 10.0;
    double inicialCondition = 1.0;
    
    if(method == 0)
    {
        ans = metodoEuler(diferential_function,ini,final,deltaT,inicialCondition);
        
    }
    else if(method == 1)
    {
        ans = rungeKutta2(diferential_function,ini,final,deltaT,inicialCondition);
    }
    else if(method == 2)
    {
        ans = rungeKutta4(diferential_function,ini,final,deltaT,inicialCondition);
    }
    else if(method == 3)
    {
        ans = taylor2Ordem(diferential_function_second,0,100,deltaT,0,0);
    }
    else if(method == 4)
    {
        ans = rungeKuttaNystron(diferential_function_second,0,20,deltaT,0,0);
    }

    std::cout << "[T]:" << std::endl;
    std::cout << "[";
    for(auto t: ans.first)
    {
        std::cout << t << ",";
        //std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "[X]:" << std::endl;
    std::cout << "[" << std::endl;

    for(auto x: ans.second)
    {
        std::cout << std::setprecision(5) << x << ",";
        //std::cout << std::endl;
    }
    std::cout << std::endl;
    divisor();
}


