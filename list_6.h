#include <vector>
#include <iostream>
#include <utility>

double derivadaPassoFrente          (double funcaoDiferencial(double x), double x, double deltaX);
double derivadaPassoTras            (double funcaoDiferencial(double x), double x,double deltaX);
double derivadaPassoCentral              (double funcaoDiferencial(double x),double x,double deltaX);
double derivadaPassoFrenteRichard   (double funcaoDiferencial(double x),double x,double deltaX,double p);
double derivadaPassoTrasRichard     (double funcaoDiferencial(double x),double x,double deltaX,double p);
double derivadaPassoCentralRichard       (double funcaoDiferencial(double x),double x,double deltaX,double p);
 