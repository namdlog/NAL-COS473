#include "Matrix.h"
#include "matrix_utils.h"
#include "Function.h"
#include <math.h>
#include <utility>
#include <vector>
#include "interface.h"
#include <algorithm>

std::pair<std::vector<double>,std::vector<double>> metodoEuler(double funcaoDiferencial(double t,double f),double tInicial,double tFinal,double deltaT,double resultadoCondicaoInicial){
    std::vector<double> incognitasX;
    incognitasX.push_back(resultadoCondicaoInicial);
    std::vector<double> incognitasT;
    incognitasT.push_back(tInicial);
    int numIteracoes = int((tFinal - tInicial)/(deltaT));
    for(int i=0;i<numIteracoes;i++){
        incognitasT.push_back((i+1)*deltaT);
        incognitasX.push_back(incognitasX[i]+deltaT*funcaoDiferencial(incognitasT[i], incognitasX[i]));
    }
    return {incognitasT,incognitasX};
}

std::pair<std::vector<double>,std::vector<double>> rungeKutta2(double funcaoDiferencial(double t,double f),double tInicial,double tFinal,double deltaT,double resultadoCondicaoInicial){
    std::vector<double> incognitasX;
    incognitasX.push_back(resultadoCondicaoInicial);
    std::vector<double> incognitasT;
    incognitasT.push_back(tInicial);
    int numIteracoes = int((tFinal - tInicial)/deltaT);
    for(int i=0;i<numIteracoes;i++){
        incognitasT.push_back((i+1)*deltaT);
        double K1 = funcaoDiferencial(incognitasT[i], incognitasX[i]);
        double K2 = funcaoDiferencial(incognitasT[i] + deltaT, incognitasX[i] + deltaT*K1);
        incognitasX.push_back(incognitasX[i] + deltaT/2*(K1 + K2));
    }
    return {incognitasT,incognitasX};
}

std::pair<std::vector<double>,std::vector<double>> rungeKutta4(double funcaoDiferencial(double t,double f),double tInicial,double tFinal,double deltaT,double resultadoCondicaoInicial){
    std::vector<double> incognitasX;
    incognitasX.push_back(resultadoCondicaoInicial);
    std::vector<double> incognitasT;
    incognitasT.push_back(tInicial);
    int numIteracoes = int((tFinal - tInicial)/deltaT);
    for(int i=0;i<numIteracoes;i++){
        incognitasT.push_back((i+1)*deltaT);
        double K1 = funcaoDiferencial(incognitasT[i], incognitasX[i]);
        double K2 = funcaoDiferencial(incognitasT[i] + deltaT/2, incognitasX[i] + deltaT/2*K1);
        double K3 = funcaoDiferencial(incognitasT[i] + deltaT/2, incognitasX[i] + deltaT/2*K2);
        double K4 = funcaoDiferencial(incognitasT[i] + deltaT, incognitasX[i] + deltaT*K3);
        incognitasX.push_back(incognitasX[i] + deltaT/6*(K1 + 2*K2+2*K3+K4));
    }
    return {incognitasT,incognitasX};
}

std::pair<std::vector<double>,std::vector<double>> taylor2Ordem(double funcaoDiferencial2Ordem(double t,double funcao, double funcaoDiferencial),double tInicial,double tFinal,double deltaT,double xInicial,double xLinhaInicial){
    std::vector<double> incognitasX;
    std::vector<double> incognitasT;
    
    incognitasX.push_back(xInicial);
    incognitasT.push_back(tInicial);
    
    double xLinhaAtual = xLinhaInicial;
    
    int numIteracoes = (int)(tFinal - tInicial)/deltaT;
    
    for(int i=0;i<numIteracoes;i++){
        incognitasT.push_back((i+1)*deltaT);
        double x2Linhas = funcaoDiferencial2Ordem(incognitasT[i], incognitasX[i], xLinhaAtual);
        incognitasX.push_back(incognitasX[i] + xLinhaAtual*deltaT + x2Linhas*deltaT*deltaT/2);
        xLinhaAtual = xLinhaAtual + x2Linhas*deltaT;
    }
  
    return {incognitasT,incognitasX};
}

std::pair<std::vector<double>,std::vector<double>> rungeKuttaNystron(double funcaoDiferencial2Ordem(double t,double funcao, double funcaoDiferencial),double tInicial,double tFinal,double deltaT,double xInicial,double xLinhaInicial){
    std::vector<double> incognitasX;
    std::vector<double> incognitasT;
    
    incognitasX.push_back(xInicial);
    incognitasT.push_back(tInicial);
    
    double xLinhaAtual = xLinhaInicial;
    
    int numIteracoes = (int)(tFinal - tInicial)/deltaT;
    
    for(int i=0;i<numIteracoes;i++){
        incognitasT.push_back((i+1)*deltaT);
        double K1          = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i], incognitasX[i], xLinhaAtual);
        double Q           = deltaT/2 * (xLinhaAtual + K1/2);
        double K2          = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i] + deltaT/2, incognitasX[i] + Q, xLinhaAtual + K1);
        double K3          = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i] + deltaT/2, incognitasX[i] + Q, xLinhaAtual + K2);
        double L           = deltaT * (xLinhaAtual + K3);
        double K4          = deltaT/2 * funcaoDiferencial2Ordem(incognitasT[i] + deltaT, incognitasX[i] + L, xLinhaAtual + 2*K3);
        incognitasX.push_back(incognitasX[i] + deltaT*(xLinhaAtual + (K1 + K2 + K3)/3));
        double xLinhaAtual = xLinhaAtual + (K1 + 2*K2 + 2*K3 + K4)/3;
    }
  
    return {incognitasT,incognitasX};
}
