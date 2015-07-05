#include <stdio.h>
#include "parallel_IO_para.h"
#include <omp.h>
int convert();
double xvalue(double x[dim1][dim2][dim3], int i, int j, int k, int p);