#pragma once
//data root path
#define root "D:\\linear\\linearSolverData\\case1\\"
#define path(exp) (root #exp)
//matrix dimention infomation
#define dim1 360
#define dim2 180
#define dim3 38
#define dim4 19
//parallel partition infomation
#define node_num 4
#define block_num1 2	//should be [sqrt(node_num)], and must be even number
#define block_num2 2	//should be node_num / block_num1
#define block_step1 (dim1/block_num1)
#define block_step2 (dim2/block_num2)
//GMRES parameters
#define k_space_dim 10	//Krylov space dimension
#define p_space_dim 10	//pre-condition matrix space dimmension
#define eps (1.0e-10)	//max solution error
#define max_iter 150	//max iteration number

//debug controler
#define post_read_data_check 1	//check 2-Norm of the initial solution to check IO correctness
#define _QUARD_FILE_ 1	//1 for use 4-file IO, 0 for 2-file IO
#define _SERIAL_SCATTER_ 0	//1 for use serial data scatter, 0 for use parallel scatter
#define _PRE_CONDITION_	 1	//1 for use pre-condition method, 0 for directly solve
#define _convert_ 0	//1 for yes, 0 for no