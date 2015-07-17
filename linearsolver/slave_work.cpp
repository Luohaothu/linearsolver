#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include "basic_operation.h"
#include "IO_functions.h"
#include "parameters.h"
#include "group.h"

void slave_work(int rank, int size, MPI_Group(&MPI_GROUP_WORLD), MPI_Group(&row_group),
	MPI_Group(&coloum_group)[block_num1], MPI_Comm(&row_comm), MPI_Comm(&coloum_comm)[block_num1],
	MPI_Datatype(&xslice), MPI_Datatype(&yslice))
{
	static double a[block_step1][block_step2][dim3][dim4];
	static double b[block_step1 + 2][block_step2 + 2][dim3];
	static double x0[block_step1 + 2][block_step2 + 2][dim3];
	static double r0[block_step1 + 2][block_step2 + 2][dim3];

#if _SERIAL_SCATTER_
	MPI_Status status[block_step1];
	MPI_Request request[block_step1];
	//receive data from Node 0
	//receive a
	for (int m = 0; m < block_step1; m++)
	{
		MPI_Irecv(&a[m][0][0][0], block_step2*dim3*dim4, MPI_DOUBLE,
			0, rank * 100 + m, MPI_COMM_WORLD, &request[m]);
	}
	int error = MPI_Waitall(block_step1, request, status);
	if (error)
	{
		printf("Error when receive A from Node 0 to Node %d : code %d\n", rank, error);
		MPI_Comm_call_errhandler(MPI_COMM_WORLD, error);
	}
	//receive x0
	for (int m = 1; m <= block_step1; m++)
	{
		MPI_Irecv(&x0[m][1][0], block_step2*dim3, MPI_DOUBLE, 0,
			rank * 100 + m, MPI_COMM_WORLD, &request[m]);
	}
	error = MPI_Waitall(block_step1, request, status);
	if (error)
	{
		printf("Error when receive x0 from Node 0 to Node %d : code %d\n", rank, error);
		MPI_Comm_call_errhandler(MPI_COMM_WORLD, error);
	}
	//receive b
	for (int m = 1; m <= block_step1; m++)
	{
		MPI_Irecv(&b[m][1][0], block_step2*dim3, MPI_DOUBLE, 0,
			rank * 100 + m, MPI_COMM_WORLD, &request[m]);
	}
	error = MPI_Waitall(block_step1, request, status);
	if (error)
	{
		printf("Error when receive b from Node 0 to Node %d : code %d\n", rank, error);
		MPI_Comm_call_errhandler(MPI_COMM_WORLD, error);
	}
#else
	//receive data from node 0
	int row_rank;
	MPI_Group_rank(row_group, &row_rank);
	MPI_Status row_status;
	MPI_Request row_request;
	MPI_Status coloum_status[block_num2 - 1];
	MPI_Request coloum_request[block_num2 - 1];
	if (row_rank != MPI_UNDEFINED_RANK)
	{
		//nodes in row group receive a from node 0
		static double a_buf[block_step1][dim2][dim3][dim4];
		static double b_buf[block_step1][dim2][dim3];
		static double x0_buf[block_step1][dim2][dim3];
		MPI_Irecv(&a_buf[0][0][0][0], block_step1*dim2*dim3*dim4,
			MPI_DOUBLE, 0, row_rank, row_comm, &row_request);
		int error = MPI_Wait(&row_request, &row_status);
		if (error)
		{
			printf("Error when receive A from Node 0, Node %d in row group members : code %d\n", row_rank, error);
			MPI_Comm_call_errhandler(row_comm, error);
		}
		//copy A to itself
		for (int i = 0; i < block_step1; i++)
		{
			for (int j = 0; j < block_step2; j++)
			{
				for (int k = 0; k < dim3; k++)
				{
					for (int m = 0; m < dim4; m++)
					{
						a[i][j][k][m] = a_buf[i][j][k][m];
					}
				}
			}
		}
		//then send a to nodes in coloum[row_rank]
		for (int m = 0; m < block_step1; m++)
		{
			for (int r = 1; r < block_num2; r++)
			{
				MPI_Isend(&a_buf[m][r*block_step2][0][0], block_step2*dim3*dim4,
					MPI_DOUBLE, r, r, coloum_comm[row_rank], &coloum_request[r - 1]);
			}
			error = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
			if (error)
			{
				printf("Error when send A to coloum group %d members : code %d\n", row_rank, error);
				MPI_Comm_call_errhandler(coloum_comm[row_rank], error);
			}
		}

		//nodes in row group receive b from node 0
		MPI_Irecv(&b_buf[0][0][0], block_step1*dim2*dim3, MPI_DOUBLE, 0, row_rank,
			row_comm, &row_request);
		error = MPI_Wait(&row_request, &row_status);
		if (error)
		{
			printf("Error when receive b from Node 0, Node %d in row group members : code %d\n", row_rank, error);
			MPI_Comm_call_errhandler(row_comm, error);
		}
		//copy b to itself
		for (int i = 0; i < block_step1; i++)
		{
			for (int j = 0; j < block_step2; j++)
			{
				for (int k = 0; k < dim3; k++)
				{
					b[i + 1][j + 1][k] = b_buf[i][j][k];
				}
			}
		}
		//then send b to nodes in coloum[row_rank]
		for (int m = 0; m < block_step1; m++)
		{
			for (int r = 1; r < block_num2; r++)
			{
				MPI_Isend(&b_buf[m][r*block_step2][0], block_step2*dim3, MPI_DOUBLE,
					r, r, coloum_comm[row_rank], &coloum_request[r - 1]);
			}
			error = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
			if (error)
			{
				printf("Error when send b to coloum group %d members : code %d\n", row_rank, error);
				MPI_Comm_call_errhandler(coloum_comm[row_rank], error);
			}
		}

		//nodes in row group receive x0 from node 0
		MPI_Irecv(&x0_buf[0][0][0], block_step1*dim2*dim3, MPI_DOUBLE, 0, row_rank,
			row_comm, &row_request);
		error = MPI_Wait(&row_request, &row_status);
		if (error)
		{
			printf("Error when recirve b from Node 0, Node %d in row group members : code %d\n", row_rank, error);
			MPI_Comm_call_errhandler(row_comm, error);
		}
		//copy x0 to itself
		for (int i = 0; i < block_step1; i++)
		{
			for (int j = 0; j < block_step2; j++)
			{
				for (int k = 0; k < dim3; k++)
				{
					x0[i + 1][j + 1][k] = b_buf[i][j][k];
				}
			}
		}
		//then send xo to nodes in coloum[row_rank]
		for (int m = 0; m < block_step1; m++)
		{
			for (int r = 1; r < block_num2; r++)
			{
				MPI_Isend(&x0_buf[m][r*block_step2][0], block_step2*dim3, MPI_DOUBLE,
					r, r, coloum_comm[row_rank], &coloum_request[r - 1]);
			}
			error = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
			if (error)
			{
				printf("Error when send x0 to coloum group %d members : code %d\n", row_rank, error);
				MPI_Comm_call_errhandler(coloum_comm[row_rank], error);
			}
		}
	}
	else
	{
		//here are sub-sub-nodes
		row_rank = rank / block_num2;
		int coloum_rank = rank % block_num2;
		//receive A from their coloum leader nodes
		for (int m = 0; m < block_step1; m++)
		{
			MPI_Irecv(&a[m][0][0][0], block_step2*dim3*dim4, MPI_DOUBLE, 0, coloum_rank,
				coloum_comm[row_rank], &coloum_request[m]);
		}
		int err = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
		if (err)
		{
			printf("Error when receive A from row leader node in coloum group %d rank %d node : code %d\n", row_rank, coloum_rank, err);
			MPI_Comm_call_errhandler(coloum_comm[row_rank], err);
		}

		//receive b from their coloum leader nodes
		for (int m = 0; m < block_step1; m++)
		{
			MPI_Irecv(&b[m + 1][0][0], block_step2*dim3, MPI_DOUBLE, 0, coloum_rank,
				coloum_comm[row_rank], &coloum_request[m]);
		}
		err = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
		if (err)
		{
			printf("Error when receive b from row leader node in coloum group %d rank %d node : code %d\n", row_rank, coloum_rank, err);
			MPI_Comm_call_errhandler(coloum_comm[row_rank], err);
		}

		//receive x0 from their coloum leader nodes
		for (int m = 0; m < block_step1; m++)
		{
			MPI_Irecv(&x0[m + 1][1][0], block_step2*dim3, MPI_DOUBLE, 0, coloum_rank,
				coloum_comm[row_rank], &coloum_request[m]);
		}
		err = MPI_Waitall(block_num2 - 1, coloum_request, coloum_status);
		if (err)
		{
			printf("Error when receive x0 from row leader node in coloum group %d rank %d node : code %d\n", row_rank, coloum_rank, err);
			MPI_Comm_call_errhandler(coloum_comm[row_rank], err);
		}
	}
#endif
	//pre-compute
	double err = 1.0;
	int iter = 0;
	int neighbor_rank[4];
	int period_neighbor_rank[2] = {MPI_PROC_NULL, MPI_PROC_NULL};
	int cart_i = rank / block_num2;
	int cart_j = rank % block_num2;
	if (cart_j == 0)
	{
		neighbor_rank[down] = MPI_PROC_NULL;
		period_neighbor_rank[down] = (rank + block_num1 / 2 * block_num2) % size;
	}
	else
	{
		neighbor_rank[down] = rank - 1;
		period_neighbor_rank[down] = MPI_PROC_NULL;
	}
	if (cart_j == block_num2 - 1)
	{
		neighbor_rank[up] = MPI_PROC_NULL;
		period_neighbor_rank[up] = (rank + block_num1 / 2 * block_num2) % size;
	}
	else
	{
		neighbor_rank[up] = rank + 1;
		period_neighbor_rank[up] = MPI_PROC_NULL;
	}
	neighbor_rank[right] = (rank + block_num2) % size;
	neighbor_rank[left] = (rank - block_num2 + size) % size;
	MPI_Barrier(MPI_COMM_WORLD);

#if _PRE_CONDITION_
	exchange_x_data(x0, rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
	//firstly, calculate the precondition polynomial
	//calculate the residual r0
	residual(a, b, x0, r0);
	double norm_r = get_global_norm(r0);
	double beta = sqrt(norm_r);
	double v[p_space_dim + 1][block_step1 + 2][block_step2 + 2][dim3] = { 0 };
	double c[p_space_dim + 1][p_space_dim + 1] = { 0 };
	double h[p_space_dim + 1][p_space_dim] = { 0 };
	sca_div_vec(r0, v[0], beta);
	c[0][0] = 1 / beta;
	//generate Hessenberg matrix
	for (int j = 0; j < p_space_dim; j++)
	{
		//first exchange v[j-1] data
		exchange_x_data(v[j], rank, neighbor_rank, period_neighbor_rank, xslice, yslice);
		//calculate init v[j]
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
	err = hessenberg_solve(h, y, beta);
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
	get_global_norm(x0);
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
		//printf("err = %2.10f, iter = %d\n", err, iter);
		iter++;
		MPI_Barrier(MPI_COMM_WORLD);
	}
}