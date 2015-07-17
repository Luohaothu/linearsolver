#pragma once
#include <mpi.h>
#include "parameters.h"
//extern MPI_Group MPI_GROUP_WORLD;
//extern MPI_Group row_group;
//extern MPI_Group coloum_group[block_num1];
//extern MPI_Comm row_comm;
//extern MPI_Comm coloum_comm[block_num1];
//extern MPI_Datatype xslice, yslice;

//neighbor rank index definition
#define up 0
#define down 1
#define left 2
#define right 3

//exchange x data function
void exchange_x_data(double(&x0)[block_step1 + 2][block_step2 + 2][dim3], int rank,
	int(&neighbor_rank)[4], int(&period_neighbor_rank)[2], MPI_Datatype(&xslice), MPI_Datatype(&yslice));

//calculate global norm function
double get_global_norm(double (&r0)[block_step1 + 2][block_step2 + 2][dim3]);

double get_global_dot_product(double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3]);

