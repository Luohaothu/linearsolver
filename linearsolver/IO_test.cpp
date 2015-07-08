#include "IO_test.h"
#define max(a,b) (a>b?a:b)
int main()
{
	//file pointer for read 
	FILE *r_fp1, *r_fp2, *r_fp3, *r_fp4;
	static double a[dim1][dim2][dim3][dim4];
	static double b[dim1][dim2][dim3];
	static double x0[dim1][dim2][dim3];
	static double p[dim1][dim2][dim3] = { 0 };
	//int err = convert();
	//if (err) return -1;
#pragma omp parallel sections
	{
#pragma omp section
		{
			r_fp1 = fopen(path(A_bi.part1), "r");
			fread(&a, sizeof(a) / 4, 1, r_fp1);
			fclose(r_fp1);
			r_fp1 = fopen(path(b_bi.part1), "r");
			fread(&b, sizeof(b) / 4, 1, r_fp1);
			fclose(r_fp1);
			r_fp1 = fopen(path(x0_bi.part1), "r");
			fread(&x0, sizeof(x0) / 4, 1, r_fp1);
			fclose(r_fp1);
		}
#pragma omp section
		{
		r_fp2 = fopen(path(A_bi.part2), "r");
		fread(&a[dim1 / 4], sizeof(a) / 4, 1, r_fp2);
		fclose(r_fp2);
		r_fp2 = fopen(path(b_bi.part2), "r");
		fread(&b[dim1 / 4], sizeof(b) / 4, 1, r_fp2);
		fclose(r_fp2);
		r_fp2 = fopen(path(x0_bi.part2), "r");
		fread(&x0[dim1 / 4], sizeof(x0) / 4, 1, r_fp2);
		fclose(r_fp2);
	}
#pragma omp section
		{
			r_fp3 = fopen(path(A_bi.part3), "r");
			fread(&a[dim1 / 2], sizeof(a) / 4, 1, r_fp3);
			fclose(r_fp3);
			r_fp3 = fopen(path(b_bi.part3), "r");
			fread(&b[dim1 / 2], sizeof(b) / 4, 1, r_fp3);
			fclose(r_fp3);
			r_fp3 = fopen(path(x0_bi.part3), "r");
			fread(&x0[dim1 / 2], sizeof(x0) / 4, 1, r_fp3);
			fclose(r_fp3);
		}
#pragma omp section
		{
			r_fp4 = fopen(path(A_bi.part4), "r");
			fread(&a[3 * dim1 / 4], sizeof(a) / 4, 1, r_fp4);
			fclose(r_fp4);
			r_fp4 = fopen(path(b_bi.part4), "r");
			fread(&b[3 * dim1 / 4], sizeof(b) / 4, 1, r_fp4);
			fclose(r_fp4);
			r_fp4 = fopen(path(x0_bi.part4), "r");
			fread(&x0[3 * dim1 / 4], sizeof(x0) / 4, 1, r_fp4);
			fclose(r_fp4);
		}
	}
	double sum = 0;
#pragma omp parallel for num_threads(8) reduction(+:sum)
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				p[i][j][k] = b[i][j][k];
				for (int m = 0; m < dim4; m++)
				{
					p[i][j][k] -= a[i][j][k][m] * xvalue(x0, i, j, k, m);
				}
				sum += p[i][j][k] * p[i][j][k];
			}
		}
	}
	printf("2-Norm : %2.10f", sum);
	scanf("%*c");
	return 0;
}