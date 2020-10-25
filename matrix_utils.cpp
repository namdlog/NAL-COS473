#include "Matrix.h"
#include "matrix_utils.h"
#include <math.h>
#include <utility>
#include "interface.h"

#define PI 3.141592653589793238462643383279502884197169399

void add(Matrix A,Matrix B,Matrix &C)
{
	for(unsigned int i = 0;i < A.get_n(); i++)
    {
        for(unsigned int j = 0;j < A.get_m(); j++)
        {
            C.entries[i][j] = A.entries[i][j] + B.entries[i][j];
        }
    }
}

void multiply(Matrix A ,Matrix B, Matrix &C)
{
    
    unsigned int n = A.get_n(), m = A.get_m(), p = B.get_n(), o = B.get_m();

    Matrix Caux(n,o,"Caux");

	for(unsigned int i = 0; i < n; ++i)
    {
		for(unsigned int j = 0; j < o; ++j)
        {
			for(unsigned int k = 0; k < m; ++k)
			{
                Caux.entries[i][j] = (Caux.entries[i][j] + A.entries[i][k] * B.entries[k][j]);
            }
        }
    }
    for(unsigned int i = 0; i < n; ++i)
    {
		for(unsigned int j = 0; j < o; ++j)
        {
            C.entries[i][j] = Caux.entries[i][j];
        }
    }

}

void multiply_scalar(double scalar, Matrix a, Matrix &b)
{
    int n = a.get_n();
    int m = a.get_m();

    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            b.entries[i][j] = scalar*a.entries[i][j];
        }
    }
}

void transpose(Matrix A,Matrix &At)
{
    unsigned int n = At.get_n();
    unsigned int m = At.get_m();
   
    for(unsigned int i = 0; i < n; i++)
    {
        for(unsigned int j = i; j < m; j++)
        {
            double aux = A.entries[i][j];
            At.entries[i][j] = A.entries[j][i];
            At.entries[j][i] = aux;
        }
    }
}

bool symmetric(Matrix A)
{
    unsigned int n = A.get_n();
    unsigned int m = A.get_m();

    for(unsigned int i = 0; i < n; i++)
    {
        for(unsigned int j = i; j < m; j++)
        {
            if(A.entries[i][j] != A.entries[j][i])
            {
                return false;
            }
        }
    }

    return true;
}

void get_cofactor(Matrix A, Matrix &temp, int p, int q, int n, int m)
{

    int i = 0, j = 0; 

    for(int row = 0; row < n; row++)
    { 
        for(int col = 0; col < m; col++)
        { 
            if(row != p && col != q)
            { 
                temp.entries[i][j++] = A.entries[row][col]; 
                if (j == n-1)
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 

bool copy_matrix(Matrix to_copy, Matrix &copy)
{
    int n = to_copy.get_n();
    int m = to_copy.get_m();
    
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            copy.entries[i][j] = to_copy.entries[i][j];
        }
    }

    return true;
}

double calculate_determinant(Matrix A)
{
    double determinant = 0.0;
    int n = A.get_n();
    int m = A.get_m();

    if (n == 1)
    {
        return A.entries[0][0]; 
    }

    Matrix auxiliary(n-1,n-1,"Determinant auxiliary"); 
  
    double sign = 1.0;
  
    for (unsigned int f = 0; f < n; f++) 
    { 
        get_cofactor(A, auxiliary, 0, f, n, m); 
        determinant += sign * A.entries[0][f] * calculate_determinant(auxiliary); 
        sign = -sign; 
    } 
  
    return determinant; 
}

void adjoint(Matrix A, Matrix &adj)
{ 
    if (adj.get_n() == 1) 
    { 
        adj.entries[0][0] = 1; 
        return; 
    } 

    double sign = 1.0;
    
    int n = A.get_n();
    int m = A.get_m();

    Matrix temporary(n-1,m-1,"Temporary"); 
  
    for (int i = 0; i < n; i++)
    { 
        for (int j = 0; j < m; j++)
        {
            get_cofactor(A, temporary, i, j, n, m); 
            sign = ((i+j)%2==0) ? 1.0: -1.0;   
            adj.entries[j][i] = (sign)*(calculate_determinant(temporary)); 
        } 
    } 
} 

bool inverse(Matrix A, Matrix &A_inverse)
{

    int n = A.get_n();
    int m = A.get_m();

    double determinant = calculate_determinant(A);

    if(determinant == 0.0) 
    { 
        std::cout << A.get_name() << " is a singular matrix, can't find its inverse" << std::endl; 
        return false;
    } 
   
    Matrix adj(n,m,"Adjoint"); 
    adjoint(A, adj); 
  
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            A_inverse.entries[i][j] = adj.entries[i][j]/double(determinant);
        }
    }
  
    return true;
}


void pivot_matrix(Matrix &Pivot, int j,int piv)
{
    Pivot.entries[j][j] = 0;
    Pivot.entries[piv][piv] = 0;
    Pivot.entries[j][piv] = 1;
    Pivot.entries[piv][j] = 1;
}

bool gauss_elimination(Matrix A, Matrix &M, Matrix &U)
{
    
    int n = A.get_n();
    int m = A.get_m();

    for(int j = 0; j < n; j++)
    {
        Matrix Maux(n,n,"M_auxiliar");
        Maux.id();
        
        if(!A.entries[j][j])
        {
            Matrix pivot(n,n,"Pivot");
            pivot.id();
            int piv = -1;
            
            for(unsigned int i = j+1; i < n; i++)
            {
                if(A.entries[i][j]){
                    piv = i;
                    break;
                }
            }
            
            if(piv == -1)
            {
                return false;
            }
            else
            {
                pivot_matrix(pivot,j,piv);
                multiply(pivot,A,A);
                multiply(pivot,M,M);
            }
        }

        for(int i = j+1; i < n; i++)
        {
            Maux.entries[i][j] = (-1)*(A.entries[i][j]/A.entries[j][j]);
        }
        
        multiply(Maux,A,A);
        multiply(Maux,M,M);

    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            U.entries[i][j] = A.entries[i][j];
        }
    }
    
    return true;
}

bool positive_definite(Matrix A){

    int n = A.get_n();
    int m = A.get_m();

    for(int i = 1; i <= n; i++)
    {
        Matrix sylvester(i,i,"Sylvester");

        for(int j=0; j < i; j++)
        {
            for(int k = 0; k < i; k++)
            {
                sylvester.entries[j][k] = A.entries[j][k];
            }
        }
        if(calculate_determinant(sylvester) <= 0)
        {
            return false;
        }
    }
        
    return true;
}

bool decomposition_lu(Matrix A, Matrix &L, Matrix &U)
{
    if(!calculate_determinant(A))
    {
        std::cout << "The decomposition in LU is not possible because the determinant is 0" << std::endl;
        return false;
    }

    int n = A.get_n();
    int m = A.get_m();

    if(n != m )
    {
        std::cout << "The decomposition in LU is not possible because the matrix is not square" << std::endl;
        return false;
    }

    Matrix M(n,m,"M");
    Matrix inv(n,m,"Inverse");
    M.id();
    inv.id();

    gauss_elimination(A,M,U);

    if(inverse(M, inv))
    {
        copy_matrix(inv,L);
    }

    return true;
}

bool decomposition_cholesky(Matrix A, Matrix &L, Matrix &Lt){
    
    int n = A.get_n();
    int m = A.get_m();

    if(!symmetric(A) && !positive_definite(A))
    {
        std::cout << "The decomposition in Cholesky is not possible because the matrix is not symetric positive-definite" << std::endl;
        return false;
    }

    if(n != m )
    {
        std::cout << "The decomposition in Cholesky is not possible because the matrix is not square" << std::endl;
        return false;
    }


    for(int i = 0; i < n; i++)
    {
        double kl = 0;
        
        for(int k = 0; k < i; k++)
        {
            kl += (L.entries[i][k] * L.entries[i][k]);
        }

        L.entries[i][i] = sqrt(A.entries[i][i]-kl);

        for(int j = i+1; j < n; j++)
        {
            kl = 0;
            
            for(int k = 0; k < i; k++)
            {
                kl += L.entries[i][k] * L.entries[j][k];
            }

            L.entries[j][i] = (1 / L.entries[i][i]) * (A.entries[i][j] - kl);
        }
    }
    
    transpose(L,Lt);

    return true;
}


void substitution_foward(Matrix A, Matrix &x, Matrix b)
{
    int n = A.get_n();

    for(int i = 0; i < n; i++)
    {
        double som = 0.0;
        for(int j = 0; j < i; j++)
        {
            som += A.entries[i][j] * x.entries[j][0];
        }
        x.entries[i][0] = (b.entries[i][0] - som) / A.entries[i][i];
    }
}

void substitution_back(Matrix A, Matrix &x, Matrix b)
{
    int n = A.get_n();

    for(int i = n-1; i >= 0; i--)
    {
        double som = 0.0;
        
        for(int j = n-1; j > i; j--)
        {
            som += A.entries[i][j] * x.entries[j][0];
        }
        x.entries[i][0] = (b.entries[i][0] - som) / A.entries[i][i];
    }
}

bool solve_linear_system(Matrix A, Matrix &x, Matrix b)
{
    int verify_substitution = 0;
    int M = A.get_m();

    for(int i=1;i<M;i++)
    {
        if(A.entries[0][i]!=0)
        {
            verify_substitution++;
        }
    }

    if(verify_substitution)
    {
        substitution_back(A,x,b);
    }else{
        substitution_foward(A,x,b);
    }
    
    return true;
}

bool solve_linear_system_trick(Matrix A, Matrix A_star, Matrix &x, Matrix b)
{
    int n = x.get_n();
    int m = x.get_m();

    Matrix y(n,m,"Y");

    solve_linear_system(A,y,b);

    solve_linear_system(A_star,x,y);
    
    return true;
}

bool diagonal_dominant(Matrix A){
    
    int n = A.get_n();

    for(int i = 0; i < n; i++)
    {

        double k = abs(A.entries[i][i]);
        double somal = 0;
        double somac = 0;
    
        for(int j = 0;j < n; j++)
        {
            if(j!=i)
            {
                somal += abs(A.entries[i][j]);
                somac += abs(A.entries[j][i]);
            }
        }
        if(k<somal or k<somac)
        {
            return false;
        }
    }
    return true;
}

double euclidean_distance(Matrix A, Matrix B)
{
    int N = A.get_n();

    double k = 0.0;

    for(int i = 0; i < N; i++)
    {
        k += abs(A.entries[i][0] - B.entries[i][0]) * abs(A.entries[i][0] - B.entries[i][0]);
    }

    return sqrt(k);
}

void jacobi(Matrix A, Matrix &x, Matrix b, double tol){
    
    int n = A.get_n();
    int m = A.get_m();
    
    Matrix y(n,1,"y");
    Matrix z(n,1,"z");

    double r = 1.0;
    
    if(!diagonal_dominant(A))
    {
        std::cout << "The matrix won't converge\n";
        return;
    }

    int step = 0;

    while(r > tol)
    {
        for(int j = 0; j < m; j++)
        {
            x.entries[j][0] = b.entries[j][0];

            for(int k = 0; k < m; k++)
            {
                if(k != j)
                {
                    x.entries[j][0] += (-1)*A.entries[j][k]*y.entries[k][0];
                }
            }
            x.entries[j][0] /= A.entries[j][j];
        }

        r = euclidean_distance(x,y)/euclidean_distance(x,z);
       
        for(int j = 0; j < m; j++)
        {
            y.entries[j][0] = x.entries[j][0];
        }
       
        // printf("Step: %d\n",step++);
        // printf("Res: %.4lf\n",r);
        // y.write();
        step++;
    
    }
    for(int i = 0; i < m; i++)
    {
        x.entries[i][0] = y.entries[i][0];
    }
    std::cout<< step << std::endl;
}

void gauss_seidel(Matrix A, Matrix &x, Matrix b,double tol){
    
    if(!diagonal_dominant(A))
    {
        if(!symmetric(A) && !positive_definite(A)){
            std::cout << "The matrix won't converge\n";
            return;
        }
    }

    
    int n = A.get_n();
    int m = A.get_m();
    
    Matrix y(n,1,"y");
    Matrix z(n,1,"z");

    double r = 1.0;
    
    int step = 0;
    
    while(r > tol)
    {
        for(int j = 0; j < m; j++)
        {
            x.entries[j][0] = b.entries[j][0];

            for(int k = 0; k < m; k++)
            {
                if(k != j)
                {
                    x.entries[j][0] += (-1) * A.entries[j][k] * x.entries[k][0];
                }
            }
            x.entries[j][0] /= A.entries[j][j];
        }

        r = euclidean_distance(x,y)/euclidean_distance(x,z);
        
        for(int j = 0; j < m; j++)
        {
            y.entries[j][0] = x.entries[j][0];
        }
        step++;
        //printf("Step: %d\n",step++);
        //printf("Res: %.4lf\n",r);
        //write_Matrix(y,"y");
    }
    std::cout << step << std::endl;
}

std::pair<int,int> indices_maior_elemento(Matrix Ak){
    double answer = 0;
    int iMax = 0,jMax = 0;
    for(int i=0;i<Ak.get_n();i++){
        for(int j=0;j<Ak.get_m();j++){
            if(i != j){
                if(abs(Ak.entries[i][j])>abs(answer)){
                    answer = Ak.entries[i][j];
                    iMax   = i;
                    jMax   = j; 
                }
            }
        }
    }
    return {iMax,jMax};
}

bool verify_jacobi(Matrix Ak,double tol)
{
    std::pair<int,int> idx = indices_maior_elemento(Ak);

    if(abs(Ak.entries[idx.first][idx.second])>tol) return true;

    return false;
}

void generate_p(Matrix Ak, Matrix &Pk){
    int n = Ak.get_n();
    
    std::pair<int,int> idxMax = indices_maior_elemento(Ak);
    
    int i = idxMax.first;
    int j = idxMax.second;
    
    double phi = 0.0;

    if(Ak.entries[i][i] != Ak.entries[j][j]){
        phi = atan(2.0 * Ak.entries[i][j] / (Ak.entries[i][i] - Ak.entries[j][j])) / 2.0;
    }else{
        phi = PI/4.0;
    }

    double sin_phi = sin(phi);  
    double cos_phi = cos(phi);
    
    Pk.entries[i][i] = cos_phi;
    Pk.entries[i][j] = -1.0*sin_phi;
    Pk.entries[j][j] = cos_phi;
    Pk.entries[j][i] = sin_phi; 

}

bool eigen_power_method(Matrix A, Matrix &xz, double &max_eigen_value, double tol){

    int n = A.get_n();
    
    max_eigen_value = 1.0;
    double laux = 1.0;

    for(int i = 0; i < n; i++)
    {
        xz.entries[i][0] = 1.0;
    }

    double r = 1;
    int cont = 0;

    while(r > tol){

        multiply(A,xz,xz);

        max_eigen_value = xz.entries[0][0];
       
        for(int i = 0; i < n; i++)
        {
            xz.entries[i][0] = xz.entries[i][0] / max_eigen_value;
        }
        
        r = abs((max_eigen_value - laux) / max_eigen_value);
        laux = max_eigen_value;
    }

    return 1;
}

bool eigen_jacobi(Matrix A,Matrix &Ak, Matrix &Xk, double tol)
{
    int n = A.get_n();
    
    if(!symmetric(A)){
        std::cout << "The Matrix is not Symetric" << std::endl;
        return 0;
    }

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            Ak.entries[i][j] = A.entries[i][j];
        }
    }

    int step = 0;
    
    while(verify_jacobi(Ak,tol)){
        Matrix Pk(n,n,"Pk");
        Pk.id();

        step++;

        generate_p(Ak,Pk);
        
        Matrix AkPk(n,n,"AkPk");
        Matrix Pkt(n,n,"Pkt");

        transpose(Pk,Pkt);
        multiply(Ak,Pk,AkPk);
        
        multiply(Pkt,AkPk,Ak);
        multiply(Xk,Pk,Xk);
        
    }
    std::cout << step << std::endl;
    return 1;
}

void mmq_calc(Matrix &x, std::vector<std::pair<double,double>> samples, std::vector<std::vector<double>> coeficients){
    int n = samples.size();
    int m = coeficients[0].size();

    Matrix P(n,m,"P");
    Matrix Pt(m,n,"Pt");
    Matrix A(m,m,"A");
    Matrix A_inv(m,m,"A_inv");
    Matrix y(n,1,"y");
    Matrix c(m,1,"c");
    
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            P.entries[i][j] = coeficients[i][j];
            Pt.entries[j][i] = coeficients[i][j];
        }
    }

    P.write();
    Pt.write();

    multiply(Pt,P,A);
    A.write();

    for(int i = 0; i < n; i++)
    {
        y.entries[i][0] = samples[i].second;
    }
    
    Pt.write();
    y.write();
    multiply(Pt,y,c);
    
    c.write();
   
    inverse(A,A_inv);
    multiply(A_inv,c,x);

}
