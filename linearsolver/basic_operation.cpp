#include <omp.h>
#include <math.h>
#include <mpi.h>
#include "IO_functions.h"
#include "parameters.h"
#include "group.h"

//calculate the dot product of (x,y)
double dot_product(double x[], double y[], long long length)
{
	double sum = 0.0;
#pragma omp parallel for num_threads(omp_get_num_procs()) reduction(+:sum)
	for (long long i = 0; i < length; i++)
	{
		sum += x[i] * y[i];
	}
	return sum;
}
//dot product (Avj, vi) version
double dot_product(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3])
{
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (int i = 1; i <= block_step1; i++)
	{
		for (int j = 1; j <= block_step2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				sum += x[i][j][k] * y[i][j][k];
			}
		}
	}
	return sum;
}
//calculate the sum vector x = x + y
void vec_add(double x[], double y[], long long length)
{
#pragma omp parallel for
	for (long long i = 0; i < length; i++)
	{
		x[i] += y[i];
	}
}

//calculate the minus vector x = x - y
void vec_minus(double x[], double y[], long long length)
{
#pragma omp parallel for
	for (long long i = 0; i < length; i++)
	{
		x[i] -= y[i];
	}
}

//calculate scaler product minus vector x = x - c * y
void sca_vec_minus(double x[], double y[], double c, long long length)
{
#pragma omp parallel for
	for (long long i = 0; i < length; i++)
	{
		x[i] -= y[i] * c;
	}
}
//sca_vec_minus A*vj - hij*vj version x = x - h * y
void sca_vec_minus(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3], double h)
{
#pragma omp parallel for
	for (int i = 0; i <= block_step1 + 1; i++)
	{
		for (int j = 0; j <= block_step2 + 1; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x[i][j][k] -= h*y[i][j][k];
			}
		}
	}
}
//calculate scaler product add vector x = x + c * y
void sca_vec_add(double x[], double y[], double c, long long length)
{
#pragma omp parallel for
	for (long long i = 0; i < length; i++)
	{
		x[i] += y[i] * c;
	}
}
//sca_vec_add x0 = x0 + V * y version
//x = x + h * y
void sca_vec_add(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3], double h)
{
	sca_vec_minus(x, y, -h);
}
//sca_div_vec x = x / h
void sca_div_vec(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&h))
{
#pragma omp parallel for
	for (int i = 0; i <= block_step1 + 1; i++)
	{
		for (int j = 0; j <= block_step2 + 1; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x[i][j][k] /= h;
			}
		}
	}
}
//sca_div_vec y = x / h
void sca_div_vec(double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3], double (&h))
{
#pragma omp parallel for
	for (int i = 0; i <= block_step1 + 1; i++)
	{
		for (int j = 0; j <= block_step2 + 1; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i][j][k] = x[i][j][k] / h;
			}
		}
	}
}

//sca_mul_vec y = x * h
void sca_mul_vec(double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3], double(&h))
{
#pragma omp parallel for
	for (int i = 0; i <= block_step1 + 1; i++)
	{
		for (int j = 0; j <= block_step2 + 1; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i][j][k] = x[i][j][k] * h;
			}
		}
	}
}

//calculate 2-norm square of a vector ||x||^2 = x .* x
double norm(double x[], long long length)
{
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (long long i = 0; i < length; i++) 
	{
		sum += x[i] * x[i];
	}
	return sum;
}
//spcified norm version for r0 square!!
double norm_r0(double(&x)[block_step1 + 2][block_step2 + 2][dim3])
{
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (int i = 1; i <= block_step1; i++)
	{
		for (int j = 1; j <= block_step2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				sum += x[i][j][k] * x[i][j][k];
			}
		}
	}
	return sum;
}

//calculate Givens factor c and s according to xi and xk
//which makes xk = 0
void givens_factor(double xi, double xk, double &c, double &s)
{
	double t = 0;
	if (xk != 0.0)
	{
		if (fabs(xk) > fabs(xi))
		{
			t = xi / xk;
			s = 1 / sqrt(1 + t*t);
			c = s * t;
		}
		else
		{
			t = xk / xi;
			c = 1 / sqrt(1 + t*t);
			s = c * t;
		}
	}
	else
	{
		c = 1.0;
		s = 0.0;
	}
}

//operate Givens transform (i, k, c, s) to a vector x
void givens_operate(double x[], long long i, long long k, const double &c, const double &s)
{
	double xi = x[i], xk = x[k];
	x[i] = c * xi + s * xk;
	x[k] = c * xk - s * xi;
}

//calculate data position in a matrix according to index
//2-dimension condition
int index(int i, int j, int ylen)
{
	return i*ylen + j;
}
//3-dimension condition
int index(int i, int j, int k, int ylen, int zlen)
{
	return (i*ylen + j)*zlen + k;
}
//4-dimension condition
int index(int i, int j, int k, int m, int jlen, int klen, int mlen)
{
	return ((i*jlen + j)*klen + k)*mlen + m;
}

//solve upper triangle matrix
void up_tri_solve(double (&R)[p_space_dim + 1][p_space_dim], double (&y)[p_space_dim], double (&b)[p_space_dim])
{
	y[p_space_dim - 1] = b[p_space_dim - 1] / R[p_space_dim - 1][p_space_dim - 1];
	for (int i = p_space_dim - 2;i >= 0;i--)
	{
		y[i] = b[i];
		for (int j = i + 1; j < p_space_dim; j++)
		{
			y[i] -= R[i][j] * y[j];
		}
		y[i] /= R[i][i];
	}
}

//solve Hessenberg linear least square equation
//H is a n+1 by n matrix
//pre-condition version
double hessenberg_solve(double(&h)[p_space_dim + 1][p_space_dim], double(&y)[p_space_dim], double(&beta))
{
	double b[p_space_dim] = { beta };
	for (int i = 1; i < p_space_dim; i++)
	{
		double c = 1.0, s = 0.0;
		givens_factor(h[i][i], h[i + 1][i], c, s);
		for (int j = 0; j < p_space_dim; j++)
		{
			double hi = h[i][j], hi1 = h[i + 1][j];
			h[i][j] = c*hi + s*hi1;
			h[i + 1][j] = c*hi1 - s*hi;
		}
		givens_operate(b, i, i + 1, c, s);
		up_tri_solve(h, y, b);
	}
	return h[p_space_dim][p_space_dim - 1];
}
//calculate residual r = b - Ax
void residual(double (&a)[block_step1][block_step2][dim3][dim4], double (&b)[block_step1 + 2][block_step2 + 2][dim3], 
	double (&x)[block_step1 + 2][block_step2 + 2][dim3], double (&r)[block_step1 + 2][block_step2 + 2][dim3])
{
#pragma omp parallel for
	for (int i = 0; i < block_step1; i++)
	{
		for (int j = 0; j < block_step2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				r[i + 1][j + 1][k + 1] = b[i + 1][j + 1][k + 1];
				for (int m = 0; m < dim4; m++)
				{
					r[i + 1][j + 1][k + 1] -= a[i][j][k][m] * xvalue(x, i, j, k, m);
					if (xvalue(x, i, j, k, m) == 0 &&
						!(k == 0 && (m == 9 || m == 10 || m == 11 || m == 12 || m == 13)) &&
						!(k == dim3 - 1 && (m == 14 || m == 15 || m == 16 || m == 17 || m == 18)))
					{
						//printf("x[%d][%d][%d][%d] == 0\n", i, j, k, m);
					}
				}
			}
		}
	}
}

//calculate y = A*x
void mat_mul_vec(double(&a)[block_step1][block_step2][dim3][dim4], double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3])
{
#pragma omp parallel for
	for (int i = 0; i < block_step1; i++)
	{
		for (int j = 0; j < block_step2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i + 1][j + 1][k] = 0;
				for (int m = 0; m < dim4; m++)
				{
					y[i + 1][j + 1][k] += a[i][j][k][m] * xvalue(x, i, j, k, m);
				}
			}
		}
	}
}
//calculate x = A*x
void mat_mul_vec(double(&a)[block_step1][block_step2][dim3][dim4], double(&x)[block_step1 + 2][block_step2 + 2][dim3])
{
	double y[block_step1 + 2][block_step2 + 2][dim3] = { 0 };
#pragma omp parallel for
	for (int i = 0; i < block_step1; i++)
	{
		for (int j = 0; j < block_step2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i + 1][j + 1][k] = 0;
				for (int m = 0; m < dim4; m++)
				{
					y[i + 1][j + 1][k] += a[i][j][k][m] * xvalue(x, i, j, k, m);
				}
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < block_step1 + 2; i++)
	{
		for (int j = 0; j < block_step2 + 2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x[i][j][k] = y[i][j][k];
			}
		}
	}
}
//calculate polynomial function of A multiply x
//which is : y = p(A)*A*x
void paax(double(&a)[block_step1][block_step2][dim3][dim4], double(&p)[p_space_dim], 
	double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3])
{
	sca_mul_vec(x, y, p[p_space_dim - 1]);
	for (int i = p_space_dim - 2; i >= 0; i--)
	{
		mat_mul_vec(a, y);
		sca_vec_add(y, x, p[i]);
	}
	mat_mul_vec(a, y);
}

//calculate b = p(A)*b
void pab(double(&a)[block_step1][block_step2][dim3][dim4], double(&p)[p_space_dim],
	double(&b)[block_step1 + 2][block_step2 + 2][dim3], int rank,
	int(&neighbor_rank)[4], int(&period_neighbor_rank)[2], MPI_Datatype(&xslice), MPI_Datatype(&yslice))
{
	double y[block_step1 + 2][block_step2 + 2][dim3] = { 0 };
	sca_mul_vec(b, y, p[p_space_dim - 1]);
	for (int i = p_space_dim - 2; i >= 0; i--)
	{
		exchange_x_data(y, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
		mat_mul_vec(a, y);
		sca_vec_add(y, b, p[i]);
	}
#pragma omp parallel for
	for (int i = 0; i < block_step1 + 2; i++)
	{
		for (int j = 0; j < block_step2 + 2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				b[i][j][k] = y[i][j][k];
			}
		}
	}
}

//calculate r = x + h * y
void xphy(double(&x)[block_step1 + 2][block_step2 + 2][dim3], double(&y)[block_step1 + 2][block_step2 + 2][dim3],
	double(&r)[block_step1 + 2][block_step2 + 2][dim3], double h)
{
#pragma omp parallel for
	for (int i = 0; i < block_step1 + 2; i++)
	{
		for (int j = 0; j < block_step2 + 2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				r[i][j][k] = x[i][j][k];
			}
		}
	}
	sca_vec_add(r, y, h);
}

//serial residual
void residual(double(&a)[dim1][dim2][dim3][dim4], double(&b)[dim1][dim2][dim3],
	double(&x0)[dim1][dim2][dim3], double(&r0)[dim1][dim2][dim3])
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				r0[i][j][k] = b[i][j][k];
				for (int m = 0; m < dim4; m++)
				{
					r0[i][j][k] -= a[i][j][k][m] * pos_x_value(x0, i, j, k, m);
				}
			}
		}
	}
}

//serial get_norm
double get_norm(double(&r)[dim1][dim2][dim3])
{
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				sum += r[i][j][k] * r[i][j][k];
			}
		}
	}
	return sum;
}

//serial sca_div_rec
void sca_div_vec(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double(&h))
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i][j][k] = x[i][j][k] / h;
			}
		}
	}
}

//serial mat_mul_vec
void mat_mul_vec(double(&a)[dim1][dim2][dim3][dim4], double(&x)[dim1][dim2][dim3],
	double(&y)[dim1][dim2][dim3])
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i][j][k] = 0;
				for (int m = 0; m < dim4; m++)
				{
					y[i][j][k] += a[i][j][k][m] * pos_x_value(x, i, j, k, m);
				}
			}
		}
	}
}

//serial dot product
double dot_product(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3])
{
	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				sum += x[i][j][k] * y[i][j][k];
			}
		}
	}
	return sum;
}

//serial sca_vec_minus x = x - h * y
void sca_vec_minus(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double h)
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x[i][j][k] -= h*y[i][j][k];
			}
		}
	}
}

//serial sca_div_vec x = x / h
void sca_div_vec(double(&x)[dim1][dim2][dim3], double h)
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x[i][j][k] /= h;
			}
		}
	}
}

//serial sca_vec_add x = x + h * y
void sca_vec_add(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double h)
{
	sca_vec_minus(x, y, -h);
}

//serial sca_mul_vec y = h * x
void sca_mul_vec(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3], double h)
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i][j][k] = x[i][j][k] * h;
			}
		}
	}
}

//serial mat_mul_vec x = A*x
void mat_mul_vec(double(&a)[dim1][dim2][dim3][dim4], double(&x)[dim1][dim2][dim3])
{
	static double y[dim1][dim2][dim3] = { 0 };
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				y[i][j][k] = 0;
				for (int m = 0; m < dim4; m++)
				{
					y[i][j][k] += a[i][j][k][m] * pos_x_value(x, i, j, k, m);
				}
			}
		}
	}
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x[i][j][k] = y[i][j][k];
			}
		}
	}
}

//serial pab b=p(A)*b
void pab(double(&a)[dim1][dim2][dim3][dim4], double(&p)[p_space_dim],
	double(&b)[dim1][dim2][dim3])
{
	static double y[dim1][dim2][dim3] = { 0 };
	sca_mul_vec(b, y, p[p_space_dim - 1]);
	for (int i = p_space_dim - 2; i >= 0; i--)
	{
		mat_mul_vec(a, y);
		sca_vec_add(y, b, p[i]);
	}
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				b[i][j][k] = y[i][j][k];
			}
		}
	}
}

//serial paax y = p(A)*A*x
void paax(double(&a)[dim1][dim2][dim3][dim4], double(&p)[p_space_dim], 
	double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3])
{
	sca_mul_vec(x, y, p[p_space_dim - 1]);
	for (int i = p_space_dim - 2; i >= 0; i--)
	{
		mat_mul_vec(a, y);
		sca_vec_add(y, x, p[i]);
	}
	mat_mul_vec(a, y);
}

//serial xphy r = x + h * y
void xphy(double(&x)[dim1][dim2][dim3], double(&y)[dim1][dim2][dim3],
	double(&r)[dim1][dim2][dim3], double h)
{
#pragma omp parallel for
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				r[i][j][k] = x[i][j][k];
			}
		}
	}
	sca_vec_add(r, y, h);
}