#pragma once
#include "parameters.h"
#include <mpi.h>
//calculate the dot product of (x,y)
double dot_product(double x[], double y[], long long length);
//dot product (Avj, vi) version
double dot_product(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3]);

//calculate the sum vector x = x + y
void vec_add(double x[], double y[], long long length);

//calculate the minus vector x = x - y
void vec_minus(double x[], double y[], long long length);

//calculate scaler product add vector x = x + c * y
void sca_vec_add(double x[], double y[], double c, long long length);
//sca_vec_add x0 = x0 + V * y version
void sca_vec_add(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3], double h);

//calculate scaler product minus vector x = x - c * y
void sca_vec_minus(double x[], double y[], double c, long long length);
//sca_vec_minus A*vj - hij*vj version x = x - h * y
void sca_vec_minus(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3], double h);
//sca_div_vec x = x / h
void sca_div_vec(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&h));
//sca_div_vec y = x / h
void sca_div_vec(double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3], double(&h));

//sca_mul_vec y = x * h
void sca_mul_vec(double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3], double(&h));

//calculate 2-norm square of a vector ||x||^2 = x .* x
double norm(double x[], long long length);
//spcified norm version for r0 square!!
double norm_r0(double(&x)[block_step1 + 2][block_step2 + 2][dim3]);

//calculate Givens factor c and s according to xi and xk
//which makes xk = 0
void givens_factor(double xi, double xk, double &c, double &s);

//operate Givens transform (i, k, c, s) to a vector x
inline void givens_operate(double x[], long long i, long long k, const double &c, const double &s);

//calculate data position in a matrix according to index
//2-dimension condition
inline int index(int i, int j, int ylen);

//3-dimension condition
inline int index(int i, int j, int k, int ylen, int zlen);

//4-dimension condition
inline int index(int i, int j, int k, int m, int jlen, int klen, int mlen);

//solve upper triangle matrix
void up_tri_solve(double(&R)[p_space_dim + 1][p_space_dim], double(&y)[p_space_dim], double(&b)[p_space_dim]);

//solve Hessenberg linear least square equation
//H is a n+1 by n matrix
double hessenberg_solve(double(&h)[p_space_dim + 1][p_space_dim], double(&y)[p_space_dim], double(&beta));

//calculate matrix multiply vector which is : x = A*x
//2-dimension condition
void mat_mul_vec(double *a, int len1, int len2, double x[]);

//calculate matrix multiply vector which is : x = A*x
//4-dimension condition
void mat_mul_vec(double *a, int len1, int len2, int len3, int len4, double x[], int lenx);

//calculate residual r = b - Ax
void residual(double(&a)[block_step1][block_step2][dim3][dim4], double(&b)[block_step1 + 2][block_step2 + 2][dim3],
	double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&r)[block_step1 + 2][block_step2 + 2][dim3]);

//calculate y = A*x
void mat_mul_vec(double(&a)[block_step1][block_step2][dim3][dim4], double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3]);
//calculate x = A*x
void mat_mul_vec(double(&a)[block_step1][block_step2][dim3][dim4], double(&x)[block_step1 + 2][block_step2 + 2][dim3]);

//calculate b = p(A)*b
void pab(double(&a)[block_step1][block_step2][dim3][dim4], double(&p)[p_space_dim],
	double(&b)[block_step1 + 2][block_step2 + 2][dim3], int rank,
	int(&neighbor_rank)[4], int(&period_neighbor_rank)[2], MPI_Datatype(&xslice), MPI_Datatype(&yslice));

//calculate polynomial function of A multiply x
//which is : y = p(A)*A*x
void paax(double(&a)[block_step1][block_step2][dim3][dim4], double(&p)[p_space_dim],
	double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3]);

//calculate r = x + h * y
void xphy(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3],
	double(&r)[block_step1 + 2][block_step2 + 2][dim3], double h);

//serial residual
void residual(double(&a)[dim1][dim2][dim3][dim4], double(&b)[dim1][dim2][dim3],
	double(&x0)[dim1][dim2][dim3], double(&r0)[dim1][dim2][dim3]);

//serial get_norm
double get_norm(double(&r)[dim1][dim2][dim3]);

//serial sca_div_rec
void sca_div_vec(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double(&h));

//serial mat_mul_vec
void mat_mul_vec(double(&a)[dim1][dim2][dim3][dim4], double(&x)[dim1][dim2][dim3],
	double(&y)[dim1][dim2][dim3]);

//serial dot product
double dot_product(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3]);

//serial sca_vec_minus x = x - h * y
void sca_vec_minus(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double h);

//serial sca_div_vec x = x / h
void sca_div_vec(double(&x)[dim1][dim2][dim3], double h);

//serial sca_vec_add x = x + h * y
void sca_vec_add(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double h);

//serial sca_mul_vec y = h * x
void sca_mul_vec(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double h);

//serial mat_mul_vec x = A*x
void mat_mul_vec(double(&a)[dim1][dim2][dim3][dim4], double(&x)[dim1][dim2][dim3]);

//serial pab b=p(A)*b
void pab(double(&a)[dim1][dim2][dim3][dim4], double(&p)[p_space_dim],
	double(&b)[dim1][dim2][dim3]);

//serial paax y = p(A)*A*x
void paax(double(&a)[dim1][dim2][dim3][dim4], double(&p)[p_space_dim],
	double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3]);

//serial xphy r = x + h * y
void xphy(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3],
	double(&r)[dim1][dim2][dim3], double h);
