#include <mpi.h>
#include "group.h"
#include "parameters.h"
void exchange_x_data(double (&x0)[block_step1 + 2][block_step2 + 2][dim3], int rank,
	int (&neighbor_rank)[4], int (&period_neighbor_rank)[2], MPI_Datatype (&xslice), MPI_Datatype(&yslice))
{
	//exchange bondary datas
	//send right yslice of data to right neighbor & receive left yslice of data from left neighbor
	MPI_Sendrecv(&x0[block_step1][0][0], 1, yslice, neighbor_rank[right], rank,
		&x0[0][0][0], 1, yslice, neighbor_rank[left], neighbor_rank[left],
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//send left yslice of data to left neighbor & receive right yslice of data from right neighbor
	MPI_Sendrecv(&x0[1][0][0], 1, yslice, neighbor_rank[left], rank,
		&x0[block_step1 + 1][0][0], 1, yslice, neighbor_rank[right], neighbor_rank[right],
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//send up xslice of data to up neighbor & receive down xslice of data from down neighbor
	MPI_Sendrecv(&x0[0][block_step2][0], 1, xslice, neighbor_rank[up], rank,
		&x0[0][0][0], 1, xslice, neighbor_rank[down], neighbor_rank[down],
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//send down xslice of data to down neighbor & receive up xslice of data from up neighbor
	MPI_Sendrecv(&x0[0][1][0], 1, xslice, neighbor_rank[down], rank,
		&x0[0][block_step2 + 1][0], 1, xslice, neighbor_rank[up], neighbor_rank[up],
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//deal with up and down period conditions separately
	if (period_neighbor_rank[up] != MPI_PROC_NULL)
	{
		MPI_Sendrecv(&x0[0][block_step2][0], 1, xslice, period_neighbor_rank[up], rank,
			&x0[0][block_step2 + 1][0], 1, xslice, period_neighbor_rank[up], period_neighbor_rank[up],
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
	if (period_neighbor_rank[down] != MPI_PROC_NULL)
	{
		MPI_Sendrecv(&x0[0][1][0], 1, xslice, period_neighbor_rank[down], rank,
			&x0[0][0][0], 1, xslice, period_neighbor_rank[down], period_neighbor_rank[down],
			MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
}