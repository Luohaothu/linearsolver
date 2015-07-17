#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "basic_operation.h"
#include "IO_functions.h"
#include "parameters.h"
#include "group.h"

void host_work(int rank, int size, MPI_Group(&MPI_GROUP_WORLD), MPI_Group(&row_group),
	MPI_Group(&coloum_group)[block_num1], MPI_Comm(&row_comm), MPI_Comm(&coloum_comm)[block_num1],
	MPI_Datatype(&xslice), MPI_Datatype(&yslice))
{
	static double a_init[dim1][dim2][dim3][dim4];
	static double b_init[dim1][dim2][dim3];
	static double x0_init[dim1][dim2][dim3];
	static double a[block_step1][block_step2][dim3][dim4];
	static double b[block_step1 + 2][block_step2 + 2][dim3];
	static double x0[block_step1 + 2][block_step2 + 2][dim3];
	static double r0[block_step1 + 2][block_step2 + 2][dim3];
	//read init data
	read_data(a_init, b_init, x0_init);
	//copy Node 0's data first
	//copy a
	for (int i = 0; i < block_step1; i++)
	{
		for (int j = 0; j < block_step2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				for (int m = 0; m < dim4; m++)
				{
					a[i][j][k][m] = a_init[i][j][k][m];
				}
			}
		}
	}
	//copy b (center part)
	for (int i = 1; i < block_step1 + 1; i++)
	{
		for (int j = 1; j < block_step2 + 1; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				b[i][j][k] = b_init[i - 1][j - 1][k];
			}
		}
	}
	//copy x0 (center part)
	for (int i = 1; i < block_step1 + 1; i++)
	{
		for (int j = 1; j < block_step2 + 1; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				x0[i][j][k] = x0_init[i - 1][j - 1][k];
			}
		}
	}
#if _SERIAL_SCATTER_
	MPI_Status status[block_step1];
	MPI_Request request[block_step1];

	//send data to other nodes
	//send a
	for (int r = 1; r < size; r++)
	{
		int i = r / block_num2;
		int j = r % block_num2;
		for (int m = 0; m < block_step1; m++)
		{
			MPI_Isend(&a_init[i*block_step1 + m][j*block_step2][0][0], block_step2*dim3*dim4, MPI_DOUBLE,
				r, r * 100 + m, MPI_COMM_WORLD, &request[m]);
		}
		int err = MPI_Waitall(block_step1, request, status);
		if (err)
		{
			printf("Error when send A from Node 0 to Node %d : code %d\n", r, err);
			MPI_Comm_call_errhandler(MPI_COMM_WORLD, err);
		}
	}
	//send x0
	for (int r = 1; r < size; r++)
	{
		int i = r / block_num2;
		int j = r % block_num2;
		for (int m = 0; m < block_step1; m++)
		{
			MPI_Isend(&x0_init[i*block_step1 + m][j*block_step2][0], block_step2*dim3, MPI_DOUBLE,
				r, r * 100 + m, MPI_COMM_WORLD, &request[m]);
		}
		int err = MPI_Waitall(block_step1, request, status);
		if (err)
		{
			printf("Error when send x0 from Node 0 to Node %d : code %d\n", r, err);
			MPI_Comm_call_errhandler(MPI_COMM_WORLD, err);
		}
	}
	//send b
	for (int r = 1; r < size; r++)
	{
		int i = r / block_num2;
		int j = r % block_num2;
		for (int m = 0; m < block_step1; m++)
		{
			MPI_Isend(&b_init[i*block_step1 + m][j*block_step2][0], block_step2*dim3, MPI_DOUBLE,
				r, r * 100 + m, MPI_COMM_WORLD, &request[m]);
		}
		int err = MPI_Waitall(block_step1, request, status);
		if (err)
		{
			printf("Error when send b from Node 0 to Node %d : code %d\n", r, err);
			MPI_Comm_call_errhandler(MPI_COMM_WORLD, err);
		}
	}
#else
	MPI_Status row_status[block_num1 - 1];
	MPI_Request row_request[block_num1 - 1];
	MPI_Status coloum_status[block_num2 - 1];
	MPI_Request coloum_request[block_num2 - 1];
	//send data to other nodes
	//send a
	//firstly, send a from node 0 to subhosts in group row_group
	for (int r = 1; r < block_num1; r++)
	{
		MPI_Isend(&a_init[r * block_step1][0][0][0], block_step1*dim2*dim3*dim4,
			MPI_DOUBLE, r, r, row_comm, &row_request[r - 1]);
	}
	int error = MPI_Waitall(block_num1 - 1, row_request, row_status);
	if (error)
	{
		printf("Error when send A to row group members : code %d\n", error);
		MPI_Comm_call_errhandler(row_comm, error);
	}
	//then send a from node 0 to nodes in group coloum_group[0]
	for (int m = 0; m < block_step1; m++)
	{
		for (int r = 1; r < block_num2; r++)
		{
			MPI_Isend(&a_init[m][r*block_step2][0][0], block_step2*dim3*dim4,
				MPI_DOUBLE, r, r, coloum_comm[0], &coloum_request[r - 1]);
		}
		error = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
		if (error)
		{
			printf("Error when send A to coloum group %d members : code %d\n", 0, error);
			MPI_Comm_call_errhandler(coloum_comm[0], error);
		}
	}

	//secondly, send b from node 0 to subhosts in group row_group
	for (int r = 1; r < block_num1; r++)
	{
		MPI_Isend(&b_init[r*block_step1][0][0], block_step1*dim2*dim3, MPI_DOUBLE,
			r, r, row_comm, &row_request[r - 1]);
	}
	error = MPI_Waitall(block_num1 - 1, row_request, row_status);
	if (error)
	{
		printf("Error when send b to row group members : code %d\n", error);
		MPI_Comm_call_errhandler(row_comm, error);
	}
	//then send b from node 0 to nodes in group coloum_group[0]
	for (int m = 0; m < block_step1; m++)
	{
		for (int r = 1; r < block_num2; r++)
		{
			MPI_Isend(&b_init[m][r*block_step2][0], block_step2*dim3, MPI_DOUBLE,
				r, r, coloum_comm[0], &coloum_request[r - 1]);
		}
		error = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
		if (error)
		{
			printf("Error when send b to coloum group %d members : code %d\n", 0, error);
			MPI_Comm_call_errhandler(coloum_comm[0], error);
		}
	}

	//lastly, send x0 from node 0 to subhosts in group row_group
	for (int r = 1; r < block_num1; r++)
	{
		MPI_Isend(&x0_init[r*block_step1][0][0], block_step1*dim2*dim3, MPI_DOUBLE,
			r, r, row_comm, &row_request[r - 1]);
	}
	error = MPI_Waitall(block_num1 - 1, row_request, row_status);
	if (error)
	{
		printf("Error when send x0 to row group members : code %d\n", error);
		MPI_Comm_call_errhandler(row_comm, error);
	}
	//then send x0 from node 0 to nodes in group coloum_group[0]
	for (int m = 0; m < block_step1; m++)
	{
		for (int r = 1; r < block_num2; r++)
		{
			MPI_Isend(&x0_init[m][r*block_step2][0], block_step2*dim3, MPI_DOUBLE,
				r, r, coloum_comm[0], &coloum_request[r - 1]);
		}
		error = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
		if (error)
		{
			printf("Error when send x0 to coloum group %d members : code %d\n", 0, error);
			MPI_Comm_call_errhandler(coloum_comm[0], error);
		}
	}

#endif
	//pre-compute
	double err = 1.0;
	int iter = 0;
	int neighbor_rank[4];
	int period_neighbor_rank[2] = {MPI_PROC_NULL, MPI_PROC_NULL};
	if (block_num2 > 1)
	{
		neighbor_rank[up] = 1;
		period_neighbor_rank[up] = MPI_PROC_NULL;
	}
	else
	{
		neighbor_rank[up] = MPI_PROC_NULL;
		period_neighbor_rank[up] = block_num1 / 2;
	}
	neighbor_rank[right] = block_num2;
	neighbor_rank[left] = (block_num1 - 1) * block_num2;
	neighbor_rank[down] = MPI_PROC_NULL;
	period_neighbor_rank[down] = (block_num1 / 2) * block_num2;
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Start computing...\n");

#if _PRE_CONDITION_
	exchange_x_data(x0, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
	//firstly, calculate the precondition polynomial
	//calculate the residual r0
	residual(a, b, x0, r0);
	double norm_r = get_global_norm(r0);
	printf("Norm : %2.10f\n", norm_r);
	double beta = sqrt(norm_r);
	double v[p_space_dim + 1][block_step1 + 2][block_step2 + 2][dim3] = { 0 };
	double c[p_space_dim + 1][p_space_dim + 1] = { 0 };
	double h[p_space_dim + 1][p_space_dim] = { 0 };
	sca_div_vec(r0, v[0], beta);
	c[0][0] = 1 / beta;
	//generate Hessenberg matrix
	for (int j = 0; j < p_space_dim; j++)
	{
		//first exchange v[j] data
		exchange_x_data(v[j], rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
		//calculate init v[j + 1]
		mat_mul_vec(a, v[j], v[j + 1]);
		for (int i = 0; i <= j; i++)
		{
			h[i][j] = get_global_dot_product(v[j + 1], v[i]);
		}
		for (int i = 0; i <= j; i++)
		{
			sca_vec_minus(v[j + 1], v[i], h[i][j]);
		}
		h[j + 1][j] = sqrt(get_global_norm(v[j + 1]));
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
	double vv[k_space_dim + 1][block_step1 + 2][block_step2 + 2][dim3] = { 0 };
	double hh[k_space_dim + 1][k_space_dim] = { 0 };
	//compute
	exchange_x_data(b, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
	//pab(a, p, b, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
	residual(a, b, x0, r0);
	norm_r = get_global_norm(r0);
	err = sqrt(norm_r);
	exchange_x_data(x0, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
	printf("norm x0: %2.10f\n", sqrt(get_global_norm(x0)));
	printf("Norm : %2.10f\n", norm_r);
	exchange_x_data(b, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
	while (err * err >= eps && iter <= max_iter)
	{
		double y[block_step1 + 2][block_step2 + 2][dim3] = { 0 };
		paax(a, p, x0, y);
		xphy(b, y, r0, -1);
		beta = sqrt(get_global_norm(r0));
		sca_div_vec(r0, vv[0], beta);
		//generate Hessenberg matrix
		for (int j = 0; j < k_space_dim; j++)
		{
			//exchange v[j] data
			exchange_x_data(vv[j], rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
			//calculate init v[j + 1]
			mat_mul_vec(a, vv[j], vv[j + 1]);
			for (int i = 0; i <= j; i++)
			{
				hh[i][j] = get_global_dot_product(vv[j + 1], vv[i]);
			}
			for (int i = 0; i <= j; i++)
			{
				sca_vec_minus(vv[j + 1], vv[i], hh[i][j]);
			}
			hh[j + 1][j] = sqrt(get_global_norm(vv[j + 1]));
			sca_div_vec(vv[j + 1], hh[j + 1][j]);
		}
		//solve Hessenberg problem
		double yy[k_space_dim] = { 0 };
		hessenberg_solve(hh, yy, beta);
		for (int i = 0; i < k_space_dim; i++)
		{
			sca_vec_add(x0, vv[i], yy[i]);
		}
		residual(a, b, x0, r0);
		norm_r = get_global_norm(r0);
		beta = sqrt(norm_r);
		err = beta;
		printf("err = %2.10f, iter = %d\n", err, iter);
		iter++;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	printf("Finished compute! error : %2.10f\n", err);
	printf("iter = %d\n", iter);
}