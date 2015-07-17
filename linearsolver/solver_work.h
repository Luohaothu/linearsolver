#pragma once

//node 0 work
void host_work(int rank, int size, MPI_Group (&MPI_GROUP_WORLD), MPI_Group (&row_group),
	MPI_Group(&coloum_group)[block_num1], MPI_Comm (&row_comm), MPI_Comm (&coloum_comm)[block_num1],
	MPI_Datatype (&xslice), MPI_Datatype (&yslice));

//other nodes work
void slave_work(int rank, int size, MPI_Group(&MPI_GROUP_WORLD), MPI_Group(&row_group),
	MPI_Group(&coloum_group)[block_num1], MPI_Comm(&row_comm), MPI_Comm(&coloum_comm)[block_num1],
	MPI_Datatype(&xslice), MPI_Datatype(&yslice));

//serial solver
void serial_solver();
