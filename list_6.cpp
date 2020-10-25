#include "Matrix.h"
#include "matrix_utils.h"
#include "Function.h"
#include <math.h>
#include <utility>
#include <vector>
#include "interface.h"
#include <algorithm>

double derivadaPassoFrente(double funcaoDiferencial(double x), double x, double deltaX)
{
    return (funcaoDiferencial(x+deltaX)-funcaoDiferencial(x))/deltaX;
}

double derivadaPassoTras(double funcaoDiferencial(double x), double x,double deltaX)
{
    return (funcaoDiferencial(x)-funcaoDiferencial(x-deltaX))/deltaX;
}

double derivadaPassoCentral(double funcaoDiferencial(double x),double x,double deltaX)
{
    return (funcaoDiferencial(x+deltaX)-funcaoDiferencial(x-deltaX))/(2*deltaX);
}

double derivadaPassoFrenteRichard(double funcaoDiferencial(double x),double x,double deltaX,double p)
{
    double d1 = derivadaPassoFrente(funcaoDiferencial,x,deltaX);
    double deltaX2 = deltaX/2;
    double d2 = derivadaPassoFrente(funcaoDiferencial,x,deltaX2);
    double q = deltaX/deltaX2;
    return d1+(d1-d2)/(pow(q,(-1*p))-1);
}

double derivadaPassoTrasRichard(double funcaoDiferencial(double x),double x,double deltaX,double p)
{
    double d1 = derivadaPassoTras(funcaoDiferencial,x,deltaX);
    double deltaX2 = deltaX/2;
    double d2 = derivadaPassoTras(funcaoDiferencial,x,deltaX2);
    double q = deltaX/deltaX2;
    return d1+(d1-d2)/(pow(q,(-1*p))-1);
}

double derivadaPassoCentralRichard(double funcaoDiferencial(double x),double x,double deltaX,double p){
    double d1 = derivadaPassoCentral(funcaoDiferencial,x,deltaX);
    double deltaX2 = deltaX/2;
    double d2 = derivadaPassoCentral(funcaoDiferencial,x,deltaX2);
    double q = deltaX/deltaX2;
    return d1+(d1-d2)/(pow(q,(-1*p))-1);
}