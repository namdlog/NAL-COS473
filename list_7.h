#include <vector>
#include <iostream>
#include <utility>

std::pair<std::vector<double>,std::vector<double>> metodoEuler  (double funcaoDiferencial           (double t,double f),double tInicial,double tFinal,double deltaT,double resultadoCondicaoInicial);
std::pair<std::vector<double>,std::vector<double>> rungeKutta2  (double funcaoDiferencial           (double t,double f),double tInicial,double tFinal,double deltaT,double resultadoCondicaoInicial);
std::pair<std::vector<double>,std::vector<double>> rungeKutta4  (double funcaoDiferencial           (double t,double f),double tInicial,double tFinal,double deltaT,double resultadoCondicaoInicial);
std::pair<std::vector<double>,std::vector<double>> taylor2Ordem (double funcaoDiferencial2Ordem     (double t,double funcao, double funcaoDiferencial),double tInicial,double tFinal,double deltaT,double xInicial,double xLinhaInicial);
std::pair<std::vector<double>,std::vector<double>> rungeKuttaNystron(double funcaoDiferencial2Ordem (double t,double funcao, double funcaoDiferencial),double tInicial,double tFinal,double deltaT,double xInicial,double xLinhaInicial);
