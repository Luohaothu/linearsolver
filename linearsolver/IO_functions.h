#include "parameters.h"

//convert ascii file to binary file, and make partition
int convert();

//return x value according to index
double pos_x_value(double x[dim1][dim2][dim3], int i, int j, int k, int p);

//read binary data (and check correctness)
int read_data(double(&a)[dim1][dim2][dim3][dim4], double(&b)[dim1][dim2][dim3], double(&x0)[dim1][dim2][dim3]);

double xvalue(double(&x)[block_step1 + 2][block_step2 + 2][dim3], int i, int j, int k, int p);