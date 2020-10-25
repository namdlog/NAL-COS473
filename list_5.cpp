#include "Matrix.h"
#include "matrix_utils.h"
#include "Function.h"
#include <math.h>
#include <utility>
#include <vector>
#include "interface.h"
#include <algorithm>

double polynomialIntegration(double function(double x),double a,double b,int nPoints)
{
    double eps  = 0.00001;
    double deltaX = (1.0*abs(b - a)) / ((1.0*nPoints - 1.0) + eps);
    std::vector<double> points;

    Matrix vetorB(nPoints,1,"Vetor B");

    for(int i = 1; i <= nPoints; i++)
    {
        points.push_back(1.0*(a+(1.0*(i-1.0))*deltaX));
        
        vetorB.entries[i-1][0] = ((pow(b,i)-pow(a,i))/i);
    }
  
    Matrix matVand(nPoints,nPoints,"Vandermonde");

    for(int i=0;i<nPoints;i++)
    {
        for(int j=0;j < nPoints; j++)
        {
            matVand.entries[i][j] = pow(points[j],i);
        }
    }

    Matrix L(nPoints,nPoints,"L");
    Matrix U(nPoints,nPoints,"U");
    decomposition_lu(matVand,L,U);
    Matrix x(nPoints,1,"x"); 
    solve_linear_system_trick(L,U,x,vetorB);
   
    double resp = 0.0;
    if(nPoints == 1)
    {
        points[0] = (a+b)/2.0;
    }
    
    for(int i = 0; i < nPoints; i++)
    {
        resp += x.entries[i][0] * function(points[i]+eps);
    }
    return resp;

}

double gaussIntegration(double function(double x),double a,double b,int nPoints,std::vector<double> W[11][2]){
    
    std::vector<double> pesos,pontos;

    double L = b-a;
    pesos    = W[nPoints][0];
    pontos   = W[nPoints][1];
    double somaIntegral = 0.0;

    for(int i=0; i < nPoints; i++)
    {
        somaIntegral += pesos[i]*function(((a+b+pontos[i]*L)/2));
    }
    return somaIntegral*(L/2);

}

/*
    
    return 0;
*/