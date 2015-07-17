#include <stdio.h>
#include <omp.h>
#include <math.h>
#include "basic_operation.h"
#include "IO_functions.h"
#include "parameters.h"

void serial_solver()
{
	static double a[dim1][dim2][dim3][dim4];
	static double b[dim1][dim2][dim3];
	static double x0[dim1][dim2][dim3];
	static double r[dim1][dim2][dim3];
	//read init data
	read_data(a, b, x0);
	//pre-compute
	double err = 1.0;
	int iter = 0;
	printf("Start computing...\n");
#if _PRE_CONDITION_
	//calculate the precondition polynomial
	residual(a, b, x0, r);
	double norm_r = get_norm(r);
	printf("Norm : %2.10f\n", norm_r);
	double beta = sqrt(norm_r);
	double v[p_space_dim + 1][dim1][dim2][dim3] = { 0 };
	double c[p_space_dim + 1][p_space_dim + 1] = { 0 };
	double h[p_space_dim + 1][p_space_dim] = { 0 };
	sca_div_vec(r, v[0], beta);
	c[0][0] = 1 / beta;
	//generate Hessenberg matrix
	for (int j = 0; j < p_space_dim; j++)
	{
		//calculate init v[j + 1]
		mat_mul_vec(a, v[j], v[j + 1]);
		for (int i = 0; i <= j; i++)
		{
			h[i][j] = dot_product(v[j + 1], v[i]);
		}
		for (int i = 0; i <= j; i++)
		{
			sca_vec_minus(v[j + 1], v[i], h[i][j]);
		}
		h[j + 1][j] = sqrt(get_norm(v[j + 1]));
		sca_div_vec(v[j + 1], h[j + 1][j]);
		//set c[x][j+1]
		c[0][j + 1] = 0;
		for (int i = 0; i <= j; i++)
		{
			c[0][j + 1] -= c[0][i] * h[i][j] / h[j + 1][j];
		}
		for (int i = 1; i <= j; i++)
		{
			c[i][j + 1] = c[i - 1][j] / h[j + 1][j];
			for (int k = 0; k <= j; k++)
			{
				c[i][j + 1] -= c[i][k] * h[k][j] / h[j + 1][j];
			}
		}
		c[j + 1][j + 1] = c[j][j] / h[j + 1][j];
	}

	//solve Hessenberg problem
	double y[p_space_dim] = { 0 };
	hessenberg_solve(h, y, beta);
	for (int i = 0; i < p_space_dim; i++)
	{
		sca_vec_add(x0, v[i], y[i]);
	}
	//form pre-condition polynomial
	double p[p_space_dim] = { 0 };
	for (int i = 0; i < p_space_dim; i++)
	{
		for (int j = i; j < p_space_dim; j++)
		{
			p[i] += c[i][j] * y[j];
		}
	}
#else
	double beta;
	double norm_r;
	double p[p_space_dim] = { 1 };
#endif
	//double beta;
	double vv[k_space_dim + 1][dim1][dim2][dim3] = { 0 };
	double hh[k_space_dim + 1][k_space_dim] = { 0 };
	//compute
	pab(a, p, b);
	residual(a, b, x0, r);
	norm_r = get_norm(r);
	err = sqrt(norm_r);
	printf("Norm x0: %2.10f\n", sqrt(get_norm(x0)));
	printf("Norm : %2.10f\n", norm_r);
	while (err * err >= eps && iter <= max_iter)
	{
		static double y[dim1][dim2][dim3] = { 0 };
		paax(a, p, x0, y);
		xphy(b, y, r, -1);
		beta = sqrt(get_norm(r));
		sca_div_vec(r, vv[0], beta);
		//generate Hessenberg matrix
		for (int j = 0; j < k_space_dim; j++)
		{
			//calculate init v[j + 1]
			mat_mul_vec(a, vv[j], vv[j + 1]);
			for (int i = 0; i <= j; i++)
			{
				hh[i][j] = dot_product(vv[j + 1], vv[i]);
			}
			for (int i = 0; i <= j; i++)
			{
				sca_vec_minus(vv[j + 1], vv[i], hh[i][j]);
			}
			hh[j + 1][j] = sqrt(get_norm(vv[j + 1]));
			sca_div_vec(vv[j + 1], hh[j + 1][j]);
		}
		//solve Hessenberg problem
		double yy[k_space_dim] = { 0 };
		hessenberg_solve(hh, yy, beta);
		for (int i = 0; i < k_space_dim; i++)
		{
			sca_vec_add(x0, vv[i], yy[i]);
		}
		residual(a, b, x0, r);
		norm_r = get_norm(r);
		beta = sqrt(norm_r);
		err = beta;
		printf("err = %2.10f, iter = %d\n", err, iter);
		iter++;
	}
	printf("Finished compute! error : %2.10f\n", err);
	printf("iter = %d\n", iter);
}
