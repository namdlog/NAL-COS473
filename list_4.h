#include <vector>
#include <iostream>
#include <utility>

double eval_function               (Function f,double x);
void   bisectionMethod             (double function(double x), double a,double b,double tol);
void   newtonMethod                (double function(double x), double x, double tol, int iterationsNumber);
void   secantMethod                (double function(double x), double xZero, double tol, int iterationsNumber);
double rootInverseInterpolation    (double function(double x), std::vector<double> x, double tol, double nInter);
void   interpolationCalculation    (std::vector<double> x, std::vector<double> y);
Matrix newtonsMethodNlq            (int n_equations, double function(int offset, vectors variables), Matrix x_zero,double tolerance,int iterationsNumber);
Matrix broydenMethodNlq            (int n_equations, double function(int offset, vectors variables), Matrix x_zero,double tolerance,int iterationsNumber);
Matrix mmqNl                       (int n_equations, double function(int offset, vectors variables), Matrix x_zero, double tolerance,int iterationsNumber);