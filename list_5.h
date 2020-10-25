#include <vector>
#include <iostream>
#include <utility>

double polynomialIntegration(double function(double x),double a,double b,int nPoints);
double gaussIntegration(double function(double x),double a,double b,int nPoints,std::vector<double> W[][2]);