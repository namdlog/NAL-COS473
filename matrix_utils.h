#include <vector>
#include <iostream>
#include <utility>

void add                        (Matrix A, Matrix B, Matrix &C);
void transpose                  (Matrix A, Matrix &At);
void multiply                   (Matrix A, Matrix B, Matrix &C);
void get_cofactor               (Matrix A, Matrix &temp, int p, int q, int n, int m);
bool symmetric                  (Matrix A);
double calculate_determinant    (Matrix A);
void adjoint                    (Matrix A, Matrix &adj);
bool inverse                    (Matrix A, Matrix &A_inverse);
bool decomposition_lu           (Matrix A,Matrix &L, Matrix &U);
bool decomposition_cholesky     (Matrix A,Matrix &L, Matrix &Lt);
bool gauss_elimination          (Matrix A, Matrix &M, Matrix &U);
void pivot_matrix               (Matrix &Pivot, int j,int piv);
bool copy_matrix                (Matrix to_copy, Matrix &copy);
void substitution_foward        (Matrix A, Matrix &x, Matrix b);
void substitution_back          (Matrix A, Matrix &x, Matrix b);
bool solve_linear_system        (Matrix A, Matrix &x, Matrix b);
double euclidean_distance       (Matrix A, Matrix B);
bool solve_linear_system_trick  (Matrix A, Matrix A_star, Matrix &x, Matrix b);
void jacobi                     (Matrix A, Matrix &x, Matrix b, double tol);
void gauss_seidel               (Matrix A, Matrix &x, Matrix b,double tol);
bool eigen_power_method         (Matrix A, Matrix &xz, double &max_eigen_value, double tol);
bool eigen_jacobi               (Matrix A,Matrix &Ak, Matrix &Xk,double tol);
void mmq_calc                   (Matrix &x, std::vector<std::pair<double,double>> samples, std::vector<std::vector<double>> coeficients);
bool positive_definite          (Matrix A);
void multiply_scalar            (double scalar, Matrix a, Matrix &b);