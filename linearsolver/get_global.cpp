#include <mpi.h>
#include "parameters.h"
#include "basic_operation.h"

double get_global_norm(double (&r0)[block_step1 + 2][block_step2 + 2][dim3])
{
	//calculate the norm of r0
	double part_norm_r = norm_r0(r0);
	double norm_r;
	MPI_Reduce(&part_norm_r, &norm_r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&norm_r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return norm_r;
}

double get_global_dot_product(double(&x)[block_step1 + 2][block_step2 + 2][dim3],
	double(&y)[block_step1 + 2][block_step2 + 2][dim3])
{
	double part = dot_product(x, y);
	double global;
	MPI_Reduce(&part, &global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return global;
}