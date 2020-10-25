#include "Matrix.h"
#include "matrix_utils.h"
#include "Function.h"
#include <math.h>
#include <utility>
#include <vector>
#include "interface.h"
#include <algorithm>

using vectors = std::vector<double>;

double derivative(double function(double x), double x)
{   
    double h = pow(10,-10);
    double top = function(x + h) - function(x);
    double bottom = h;
    double answer = top/bottom;

    return answer;
}

double parcialDerivative(double function(int offset, vectors variables), int offset, Matrix x_zero, int relation)
{
    double h        = pow(10,-10);
    
    int n           = x_zero.get_n();
    vectors var(n);

    for(int i=0;i<n;i++)
    {
        var[i] = x_zero.entries[i][0];
    }

    double aux      = function(offset, var);    

    var[relation]  += h;
    double top      = function(offset, var) - aux;
    double bottom   = h;
    return top/bottom;
}

void bisectionMethod(double function(double x), double a,double b, double tol)
{
    double x,fi;
    double fa,fb;

    fa = function(a);

    fb = function(b);

    if(fa < 0 && fb > 0 or fa > 0 && fb < 0)
    {
        //continue;    
    }else{
        std::cout << "Cant find bissection method root for this limits." << std::endl;
        return;
    }
    int count = 0;
    while(abs(b - a) > tol)
    {
        x = (a+b) / 2.0;
        
        std::cout << count++ << " x = " << x << std::endl; 

        fi = function(x);
        if(fi > 0.0)
        {
            b = x;
        }else
        {
            a = x;
        }
    }
    
    std::cout << "The root is: " << x << std::endl;
    return;
}

void newtonMethod(double function(double x), double x, double tol, int iterationsNumber){
    double tolk,xAux;
    double fx,flx;
    for(int i = 0; i < iterationsNumber; i++)
    {
        std::cout << i << " x = " << x << std::endl;
        xAux = x;
        fx   = function(xAux);
        flx  = derivative(function,xAux);
        x    = x - (fx / flx);
        tolk = abs(x-xAux);
        
        if(tolk < tol)
        {
            std::cout << "The root is " << x << std::endl;
            return;
        }
    }

    std::cout << "Cant converge" << std::endl;
    return;
}

void secantMethod(double function(double x), double xZero, double tol, int iterationsNumber){
    double tolk;
    double deltaX = 0.001;
    double xAnterior = xZero;
    double xAtual = xAnterior + deltaX;
    double fa = function(xZero);

    double fi;
    double xProximo;
    for(int i = 0; i < iterationsNumber; i++)
    {
        fi       = function(xAtual);
        xProximo = xAtual - fi*(xAtual-xAnterior)/(fi-fa);
        tolk     = abs(xProximo - xAtual);
        std::cout << i << " x = " << xProximo << std::endl;
        if(tolk < tol)
        {
            std::cout << "The root is " << xProximo << std::endl;
            return;
        
        }else
        {
            xAnterior = xAtual;
            xAtual = xProximo;
            fa = fi;
        }
    }

    std::cout << "Cant converge" << std::endl;
    return;
}

double interpolationCalculation(std::vector<double> x, std::vector<double> y){
    double answer = 0.0;
    for(int i=0;i<x.size();i++)
    {
        double num = 1.0;
        double den = 1.0;
        for(int j=0;j<x.size();j++){
            if(j!=i){
                num *= y[j];
            }   
        }
        for(int j=0;j<x.size();j++){
            if(j!=i){
                den *= (y[i]-y[j]);
            }    
        }
        answer += (num/den)*x[i];
    }
    return answer;
}

void rootInverseInterpolation(double function(double x), std::vector<double> x, double tol, double iterationsNumber){
    
    double xAtual,xAnterior = pow(10,36);
    int nSamples = x.size();
    std::vector<double> y(nSamples);

    for(int k = 1;k <= iterationsNumber; k++)
    {
        for(int i = 0; i < x.size(); i++)
        {
            y[i] = function(x[i]);
            //std::cout << "x -> " << x[i] << std::endl;
            //std::cout << "y -> " << y[i] << std::endl;
        }
        
        xAtual = interpolationCalculation(x,y);
        std::cout << k << " x* : " << xAtual << std::endl;
        double tolk = abs(xAtual-xAnterior);
        
        if(tolk < tol)
        {
            std::cout << "The root is: " << xAtual << std::endl;
            return;
        }else
        {
            int idxMax;
            double yMax = 0;
            for(int i=0;i<y.size();i++){
                if(abs(y[i])>yMax){
                    idxMax = i;
                    yMax = abs(y[i]);
                }
            }
            
            x[idxMax] = xAtual;
            sort(x.begin(),x.end());
            y[idxMax] = function(xAtual);
            sort(y.begin(),y.end());
            std::cout << "x* -> " << xAtual << std::endl;
            std::cout << "tolk -> " << tolk << std::endl;
            xAnterior = xAtual;
        }
    }

    std::cout << "Cant converge." << std::endl;
    return;
}

Matrix calculateJ(int n, int m, double function(int offset, vectors variables), Matrix x)
{
    Matrix J(n,m,"J");
    //J.write();
    std::vector<double> var(m,0);

    for(int i=0;i<m;i++)
    {
        var[i] = x.entries[i][0];
    }
   
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            double k = parcialDerivative(function,i,x,j);
            J.entries[i][j] = k;
        }
    }
    return J;
}

Matrix calculateF(int n, double function(int offset, vectors variables), Matrix x)
{
    Matrix F(n,1,"F");
    
    std::vector<double> var(n,0);

    for(int i=0;i<n;i++)
    {
        var[i] = x.entries[i][0];
    }
    
    for(int i=0;i<n;i++){
        std::cout << "f: " << var[i] << " ";
    }
    std::cout << std::endl;

    for(int i = 0; i < n ; i++)
    {
        F.entries[i][0] = function(i,var);
    }

    F.write();
    return F;
}

Matrix newtonsMethodNlq(int n_equations, double function(int offset, vectors variables), Matrix x_zero,double tolerance,int iterationsNumber)
{
    int n = n_equations;
    int m = n_equations;

    Matrix x(n,1,"x_solution");
    copy_matrix(x_zero,x);
    Matrix x_delta(n,1,"x_delta");
    Matrix zero(n,1,"Zero");
    Matrix F(n,1,"Function");
    Matrix J(n,m,"Jacobian");
    Matrix J_inverse(n,m,"Jacobian Inverse");
    Matrix note(1,1,"note");
    note.entries[0][0] = -1;
    
    double tolk;

    for(int i = 0; i < iterationsNumber; i++)
    {
        std::cout << i << " " << std::endl;

        J         = calculateJ(n,m,function,x);
        F         = calculateF(n,function,x);
        
        J.write();
        F.write();

        inverse(J,J_inverse);
        multiply(F,note,F); 
        multiply(J_inverse,F,x_delta);
        add(x,x_delta,x);
        
        x.write();
        
        tolk      = euclidean_distance(x_delta,zero)/euclidean_distance(x,zero);
        
        //std::cout << tolk << std::endl;
        
        if(tolk < tolerance)
        {
            return x;
        }
    }
    
    std::cout << "Convergence not reached" << std::endl;
    return x;
}

Matrix division(Matrix A,Matrix B)
{
    int nn = A.get_n();
    int n  = B.get_n();
    int m  = B.get_m();
    Matrix bInv(n,m,"bInverse");
    inverse(B,bInv);
    Matrix resp(nn,m,"Answer");
    multiply(A,bInv,resp);
    return resp;
}

Matrix broydenMethodNlq(int n_equations, double function(int offset, vectors variables), Matrix x_zero,double tolerance,int iterationsNumber)
{
    int n = n_equations;
    int m = n_equations;

    Matrix x(n,1,"x_solution");
    copy_matrix(x_zero,x);
    Matrix x_delta(n,1,"x_delta");
    Matrix zero(n,1,"Zero");
    Matrix F_old_x(n,1,"Old function");
    Matrix F_new_x(n,1,"New function");
    Matrix y(n,1,"Y");
    Matrix F(n,1,"Function");
    Matrix J(n,m,"Jacobian");
    Matrix J_inverse(n,m,"Jacobian Inverse");
    Matrix B(n,m,"B");
    Matrix note(1,1,"note");
    Matrix id_negative(n,m,"id_negative");
    
    for(int i=0;i<n;i++){for(int j=0;j<m;j++){if(i == j){id_negative.entries[i][j] = -1;}else{id_negative.entries[i][j] = 0;}}}
    
    note.entries[0][0] = -1;
    B                  = calculateJ(n,m,function,x);
    
    double tolk;

    for(int i = 0; i < iterationsNumber; i++)
    {
        F_old_x = calculateF(n,function,x);
        F       = F_old_x;
        J       = B;
    
        J.write();
        F.write();

        inverse(J,J_inverse);
        multiply(F,note,F); 
        multiply(J_inverse,F,x_delta);
        
        add(x,x_delta,x);
        
        F_new_x = calculateF(n,function,x);

        multiply(F_old_x,note,F_old_x);
        add(F_new_x,F_old_x,y);
        
        x.write();
        y.write();

        tolk      = euclidean_distance(x_delta,zero)/euclidean_distance(x,zero);
        std::cout << tolk << std::endl;
        if(tolk < tolerance)
        {
            return x;
        }
        else
        {
            B.write();
            y.write();
            x_delta.write();
            //x_delta_t.write();

            Matrix x_delta_t(x_delta.get_m(),x_delta.get_n(),"X Delta Transpose");
            transpose(x_delta,x_delta_t);
            Matrix T(B.get_n(),B.get_m(),"T");
            Matrix t(B.get_n(),x_delta.get_m(),"pre-t");
            Matrix bot(1,1,"bottom");
            
            multiply(x_delta_t,x_delta,bot);
            //bot.write();
            double Bottom = bot.entries[0][0];
            multiply(B,x_delta,t);
            //t.write();
            multiply_scalar(-1,t,t);
            //y.write();
            //t.write();
            add(y,t,t);
            //t.write();
            //x_delta_t.write();
            multiply(t,x_delta_t,T);
            //T.write();
            multiply_scalar(1/Bottom,T,T);
            add(B,T,B);
            //B.write();
            //y.write();
            //x_delta.write();
            //x_delta_t.write();
        }
    }
    
    std::cout << "Convergence not reached" << std::endl;
    return x;
}

Matrix mmqNl(int n_equations, double function(int offset, vectors variables), Matrix x_zero,double tolerance,int iterationsNumber)
{
    
    int n = n_equations;
    int m = x_zero.get_n();
    Matrix B(n,1,"B");
    Matrix delta_B(n,1,"Delta B");
    Matrix J(m,n,"Jacobian");
    Matrix Jt(n,m,"Jacobian Transpose");
    Matrix F(m,1,"F");
    copy_matrix(x_zero,B);

    Matrix zero(n,1,"zero");

    for(int i = 0; i < iterationsNumber; i++)
    {
        B.write();
        J = calculateJ(n,m,function,B);
        J.write();
        transpose(J,Jt);
        Jt.write();
        F = calculateF(n,function,B);
        F.write();
        Matrix a(J.get_n(),J.get_m(),"a");
        a.write();
        multiply(Jt,J,a);
        a.write();
        inverse(a,a);
        J.write();
        a.write();
        /*
        if(a == 0)
        {
            std::cout << "Erro" << std::endl;
        }
        */
        Matrix b(a.get_n(),a.get_m(),"b");
        Matrix c(J.get_m(),1,"c");
    
        multiply_scalar(-1,a,b);
        b.write();
        multiply(Jt,F,c);
        c.write();
        multiply(b,c,delta_B);
        delta_B.write();
        add(B,delta_B,B);
        B.write();
        double tolk = euclidean_distance(delta_B,zero)/euclidean_distance(B,zero);
        if(tolk < tolerance)
        {
            return B;
        }
    }
    std::cout << "Convergence not reached" << std::endl;
    return x_zero;
}

/*
double eval_function(Function f,double x){

    
    double resp = 0.0;
    for(int i=0;i<=f.grau;i++){
        
        double val = f.coef[i]*pow(x,i);
        if(f.cos[i]){
            val *= cos(x);
        }
        resp += val;
     
    }
    
    return resp;
}

double bisectionMethod(Function f,double a,double b,double tol)
{
    double x,fi;

    while(abs(b - a) > tol)
    {
        x = (a+b) / 2.0;
        fi = eval(f,x);
        if(fi > 0.0)
        {
            b = x;
        }else
        {
            a = x;
        }
    }
    return x;
}

double raizNewton(func f, func fl, double x, double tol, int nInter){
    double tolk,xaux;
    for(int i=0;i<nInter;i++){
        xaux = x;
        x = x-(eval(f,x)/eval(fl,x));
        tolk = abs(x-xaux);
        if(tolk<tol){
            return x;
        }
    }
    printf("Não converge");
    return 0.0;
}

double interpolationCalculation(vector<double> x, vector<double> y){
    double a = (y[1]*y[2]*x[0])/((y[0]-y[1])*(y[0]-y[2]));
    double b = (y[0]*y[2]*x[1])/((y[1]-y[0])*(y[1]-y[2]));
    double c = (y[0]*y[1]*x[2])/((y[2]-y[0])*(y[2]-y[1]));
    return a + b + c;
}

double raizInterpolInv(func f, vector<double> x, double tol, double nInter){
    
    double xAtual,xAnterior = pow(10,36);
    vector<double> y(3);

    for(int k=1;k<=nInter;k++){
        for(int i=0;i<x.size();i++){
            y[i] = eval(f,x[i]);
            cout << "x -> " << x[i] << endl;
            cout << "y -> " << y[i] << endl;
        }
        
        xAtual = interpolationCalc(x,y);
        cout << xAtual << endl;
        double tolk = abs(xAtual-xAnterior);
        if(tolk < tol){
            return xAtual;
        }else{
            int idxMax;
            double yMax = 0;
            for(int i=0;i<y.size();i++){
                if(abs(y[i])>yMax){
                    idxMax = i;
                    yMax = abs(y[i]);
                }
            }
            
            x[idxMax] = xAtual;
            sort(x.begin(),x.end());
            y[idxMax] = eval(f,xAtual);
            sort(y.begin(),y.end());
            cout << "x* -> " << xAtual << endl;
            cout << "tolk -> " << tolk << endl;
            xAnterior = xAtual;
        }
    }
    printf("Convergencia não alcançada");
    return 0.0;
}

mat J(vector<funcVarVar> fs,mat b){
    int n = fs.size();
    int m = b.size();
    mat ret = zeros(n,m);
    for(int i=0;i<fs.size();i++){
        for(int j=0;j<fs[i].vars.size();j++){
            func derivadaParcial;
            zeraFunc(derivadaParcial);
            func derivar = fs[i].vars[j];
            derivadaParcial.grau = derivar.grau;
            for(int k = 0;k < derivadaParcial.grau; k++){
                derivadaParcial.coef[k] = (k+1)*derivar.coef[k+1];
            }
            double fg = eval(derivadaParcial,b[j][0]);
            ret[i][j] = fg;
        }
    }
    
    return ret;
}

mat F(vector<funcVarVar> fs,mat b){
    int n = fs.size();
    int m = b.size();
    mat ret = zeros(n,1);
    for(int i=0;i<n;i++){
        double val = 0.0;
        for(int j=0;j<m;j++){
            val += eval(fs[i].vars[j],b[j][0]);
        }
        ret[i][0] = val;
    }
    return ret;
}

mat add(mat a,mat b){
    int n = a.size();
    int m = b[0].size();
    mat resp = zeros(n,m);
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            resp[i][j] = a[i][j] + b[i][j];
        }
    }
    return resp;
}

double norma(mat a){
    double resp = 0.0;
    int n = a.size();
    for(int i=0;i<n;i++){
        resp += a[i][0]*a[i][0];
    }
    resp = sqrt(resp);
    return resp;
}

mat metodoNewtonEqNl(vector<funcVarVar> fs,mat xAtual,double tol,int nInter){
    int n = fs.size();
    int m = xAtual.size();
    
    mat nossaX = zeros(m,1);
    for(int i=0;i<m;i++){
        nossaX[i][0] = xAtual[i][0];
    }

    double tolk = 0.0;
    for(int k=0;k<nInter;k++){
        mat inversJacob = zeros(n,m);
        mat jacobMat = zeros(n,m);
        mat vF = zeros(n,1);
        mat deltaX = zeros(m,1);
       
  
        jacobMat = J(fs,nossaX);
        vF = F(fs,nossaX);
        inverse(jacobMat,inversJacob);
        inversJacob = mul_escalar(-1.0,inversJacob);
        
        deltaX = mul(inversJacob,vF);
   

       
        nossaX = add(nossaX,deltaX);
        
        double a = norma(nossaX);
        double b = norma(deltaX);
            
        tolk = b/a;
       
        if(tolk < tol){
            return nossaX;
        }
    }
    printf("Convergência não alcançada");
    return nossaX;
}

mat div(mat a,mat b){
    int nn = a.size();
    int n = b.size();
    int m = b[0].size();
    mat bInv = zeros(n,m);
    inverse(b,bInv);
    mat resp = zeros(nn,m);
    resp = mul(a,bInv);
    return resp;
}

mat B(mat jac,mat y,mat del){
    write_Matrix(jac,"B- JAC");
    write_Matrix(y,"B- Y");
    write_Matrix(del,"B- DEL");
    mat ba = mul(jac,del);
    ba = mul_escalar(-1.0,ba);
    mat be = add(y,ba);
    mat delt = trans(del);
    mat sime = mul(delt,del);
    mat bi = mul(be,delt);
    mat bo = mul_escalar(1.0/sime[0][0],bi);
    write_Matrix(bo,"BO");
    mat resp = add(jac,bo);
    write_Matrix(resp,"RESP");
    return resp;
}

mat metodoBroydenEqNl(vector<funcVarVar> fs,mat xAtual,double tol,int nInter){
    int n = fs.size();
    int m = xAtual.size();
    
    mat nossaX = zeros(m,1);
    for(int i=0;i<m;i++){
        nossaX[i][0] = xAtual[i][0];
    }
    
    mat jacobMat = J(fs,nossaX);

    double tolk = 0.0;
    for(int k=0;k<nInter;k++){
        mat inversJacob = zeros(n,m);
        mat vF = zeros(n,1);
        mat deltaX = zeros(m,1);
        

        vF = F(fs,nossaX);
       
        inverse(jacobMat,inversJacob);
        inversJacob = mul_escalar(-1.0,inversJacob);
        
        deltaX = mul(inversJacob,vF);
       
        nossaX = add(nossaX,deltaX);
        
        double a = norma(nossaX);
        double b = norma(deltaX);
            
        tolk = b/a;
       
        if(tolk < tol){
            return nossaX;
        }

        mat vFaux = zeros(n,1);
        mat Y = zeros(n,1);
        vFaux = F(fs,nossaX);
        vF = mul_escalar(-1.0,vF);
        Y = add(vF,vFaux);
        jacobMat = B(jacobMat,Y,deltaX);
    }
    printf("Convergência não alcançada");
    return nossaX;
}*/