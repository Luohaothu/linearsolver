#include <stdio.h>
#include <mpi.h>
#include "parameters.h"
#include "solver_work.h"

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	//for serial condition, use standalone function to solve
	//else start scatter initial data
	if (size == 1)
	{
		serial_solver();
		MPI_Finalize();
		return 0;
	}
	else 
	{
		MPI_Group MPI_GROUP_WORLD;
		MPI_Group row_group;
		MPI_Group coloum_group[block_num1];
		MPI_Comm row_comm;
		MPI_Comm coloum_comm[block_num1];
#if !_SERIAL_SCATTER_
		//get overall group
		MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);

		//creat row group

		//creat coloum group array

		//set row group
		int row_range[1][3] = { { 0, size - 1, block_num2 } };
		MPI_Group_range_incl(MPI_GROUP_WORLD, 1, row_range, &row_group);

		//set coloum group array
		for (int i = 0; i < block_num1; i++)
		{
			int coloum_range[1][3] = { { i*block_num2, (i + 1)*block_num2 - 1, 1 } };
			MPI_Group_range_incl(MPI_GROUP_WORLD, 1, coloum_range, &coloum_group[i]);
		}

		//creat row communicator
		MPI_Comm_create(MPI_COMM_WORLD, row_group, &row_comm);

		//creat coloum communicator array
		for (int i = 0; i < block_num1; i++)
		{
			MPI_Comm_create(MPI_COMM_WORLD, coloum_group[i], &coloum_comm[i]);
		}
#endif
		MPI_Datatype xslice, yslice;
		MPI_Type_vector(block_step1 + 2, dim3, (block_step2 + 2) * dim3, MPI_DOUBLE, &xslice);
		MPI_Type_vector((block_step2 + 2) * dim3, 1, 1, MPI_DOUBLE, &yslice);
		MPI_Type_commit(&xslice);
		MPI_Type_commit(&yslice);
		if (rank == 0)
		{
			host_work(rank, size, MPI_GROUP_WORLD, row_group,
				coloum_group, row_comm, coloum_comm,
				xslice, yslice);
		}
		else
		{
			slave_work(rank, size, MPI_GROUP_WORLD, row_group,
				coloum_group, row_comm, coloum_comm,
				xslice, yslice);
		}
		MPI_Finalize();
		return 0;
	}
}