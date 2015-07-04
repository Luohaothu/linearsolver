#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "parallel_IO_para.h"

#define max(a,b) ((a>b)?a:b)
#define abs(a) (((a)>0)?(a):-(a))
int main()
{
	static double a[dim1][dim2][dim3][dim4];
	static double b[dim1][dim2][dim3][dim4];
	//file pointer for read ascii data
	FILE *r_fp;
	//file pointer for write ascii data
	FILE *aw_fp1, *aw_fp2, *aw_fp3, *aw_fp4;
	//file pointer for read binary data
	FILE *r_fp1, *r_fp2, *r_fp3, *r_fp4;
	//file pointer for write binary data
	FILE *w_fp1, *w_fp2, *w_fp3, *w_fp4;
	clock_t t1 = clock();
	if ((r_fp = fopen(path(A.txt), "r")) == NULL)
	{
		puts("Failed open A.txt");
		return 0;
	}
	for (int i = 0; i < dim1 * dim2 * dim3 * dim4; i++)
	{
		short int w, x, y, z;
		fscanf(r_fp, "%d%d%d%d", &w, &x, &y, &z);
		if (w >= dim1 || w < 0 || x >= dim2 || x < 0 || y >= dim3 || y < 0 || z >= dim4 || z < 0)
		{
			puts("Error bound\n");
			return 0;
		}
		else
		{
			fscanf(r_fp, "%lf\n", &(a[w][x][y][z]));
		}
		if (i % (dim2 * dim3 * dim4) == 0)
		{
			printf("\r%5.2f%% Finished read", (double)i / (dim1 * dim2 * dim3 * dim4) * 100);
		}
	}
	printf("\rFinished total reading.\n");
	clock_t t2 = clock();
	printf("Time for reading A : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	fclose(r_fp);
#if _QUARD_FILE_
	if (	(w_fp1 = fopen(path(A_bi.part1), "wb")) == NULL
		||	(w_fp2 = fopen(path(A_bi.part2), "wb")) == NULL
		||  (w_fp3 = fopen(path(A_bi.part3), "wb")) == NULL
		||  (w_fp4 = fopen(path(A_bi.part4), "wb")) == NULL)
	{
		puts("Failed open A_bi");
		return 0;
	}
	fwrite(&a, sizeof(a) / 4, 1, w_fp1);
	fwrite(&(a[dim1 / 4]), sizeof(a) / 4, 1, w_fp2);
	fwrite(&(a[dim1 / 2]), sizeof(a) / 4, 1, w_fp3);
	fwrite(&(a[3 * dim1 / 4]), sizeof(a) / 4, 1, w_fp4);
	printf("Successfully convert to binary file\n");
	t1 = clock();
	printf("Time for write : %f s\n", (double)(t1 - t2) / CLOCKS_PER_SEC);
	fclose(w_fp1);
	fclose(w_fp2);
	fclose(w_fp3);
	fclose(w_fp4);
	t1 = clock();
	if (   (r_fp1 = fopen(path(A_bi.part1), "r")) == NULL
		|| (r_fp2 = fopen(path(A_bi.part2), "r")) == NULL
		|| (r_fp3 = fopen(path(A_bi.part3), "r")) == NULL
		|| (r_fp4 = fopen(path(A_bi.part4), "r")) == NULL)
	{
		puts("Failed open A_bi");
		return 0;
	}
	fread(&b, sizeof(b) / 4, 1, r_fp1);
	fread(&b[dim1 / 4], sizeof(b) / 4, 1, r_fp2);
	fread(&b[dim1 / 2], sizeof(b) / 4, 1, r_fp3);
	fread(&b[3 * dim1 / 4], sizeof(b) / 4, 1, r_fp4);
	fclose(r_fp1);
	fclose(r_fp2);
	fclose(r_fp3);
	fclose(r_fp4);
	t2 = clock();
	printf("Time for read A_bi : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	if ((aw_fp1 = fopen(path(A_ascii.part1), "w")) == NULL
		|| (aw_fp2 = fopen(path(A_ascii.part2), "w")) == NULL
		|| (aw_fp3 = fopen(path(A_ascii.part3), "w")) == NULL
		|| (aw_fp4 = fopen(path(A_ascii.part4), "w")) == NULL)
	{
		puts("Failed open A_ascii");
		return 0;
	}
	t2 = clock();
#pragma omp sections
	{
#pragma omp section		//write thread 1
	for(int i = 0; i < dim1 / 4; i++)
	for(int j = 0; j < dim2; j++)
	for(int k = 0; k < dim3; k++)
	for(int m = 0; m < dim4; m++)
		fprintf(aw_fp1, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);

#pragma omp section		//write thread 2
	for(int i = dim1 / 4; i < dim1 / 2; i++)
	for(int j = 0; j < dim2; j++)
	for(int k = 0; k < dim3; k++)
	for(int m = 0; m < dim4; m++)
		fprintf(aw_fp2, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);

#pragma omp section		//write thread 3
	for(int i = dim1 / 2; i < 3 * dim1 / 4; i++)
	for(int j = 0; j < dim2; j++)
	for (int k = 0; k < dim3; k++)
	for (int m = 0; m < dim4; m++)
		fprintf(aw_fp3, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);

#pragma omp section		//write thread 4
	for (int i = 3 * dim1 / 4; i < dim1; i++)
	for (int j = 0; j < dim2; j++)
	for (int k = 0; k < dim3; k++)
	for (int m = 0; m < dim4; m++)
		fprintf(aw_fp4, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);
	}
	t1 = clock();
	printf("Time for write ascii file: %f s\n", (double)(t1 - t2) / CLOCKS_PER_SEC);
	fclose(aw_fp1);
	fclose(aw_fp2);
	fclose(aw_fp3);
	fclose(aw_fp4);
#else
	if ((w_fp1 = fopen(path(A_bi.part1), "wb")) == NULL
		|| (w_fp2 = fopen(path(A_bi.part2), "wb")) == NULL)
	{
		puts("Failed open A_bi");
		return 0;
	}
	fwrite(&a, sizeof(a) / 2, 1, w_fp1);
	fwrite(&(a[dim1 / 2]), sizeof(a) / 2, 1, w_fp2);
	printf("Successfully convert to binary file\n");
	t1 = clock();
	printf("Time for write : %f s", (double)(t1 - t2) / CLOCKS_PER_SEC);
	fclose(w_fp1);
	fclose(w_fp2);
	t1 = clock();
	if ((r_fp1 = fopen(path(A_bi.part1), "r")) == NULL
		|| (r_fp2 = fopen(path(A_bi.part2), "r")) == NULL)
	{
		puts("Failed open A_bi");
		return 0;
	}
	fread(&b, sizeof(b) / 2, 1, r_fp1);
	fread(&b[dim1 / 2], sizeof(b) / 2, 1, r_fp2);
	fclose(r_fp1);
	fclose(r_fp2);
	t2 = clock();
	printf("Time for read A_bi : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	if ((aw_fp1 = fopen(path(A_ascii.part1), "w")) == NULL
		|| (aw_fp2 = fopen(path(A_ascii.part2), "w")) == NULL)
	{
		puts("Failed open A_ascii");
		return 0;
	}
	t2 = clock();
#pragma omp sections
	{
#pragma omp section		//write thread 1
		for (int i = 0; i < dim1 / 2; i++)
		for (int j = 0; j < dim2; j++)
		for (int k = 0; k < dim3; k++)
		for (int m = 0; m < dim4; m++)
			fprintf(aw_fp1, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);

#pragma omp section		//write thread 2
		for (int i = dim1 / 2; i < dim1; i++)
		for (int j = 0; j < dim2; j++)
		for (int k = 0; k < dim3; k++)
		for (int m = 0; m < dim4; m++)
			fprintf(aw_fp2, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);
	}
	fclose(aw_fp1);
	fclose(aw_fp2);
	t1 = clock();
	printf("Time for write ascii file: %f s\n", (double)(t1 - t2) / CLOCKS_PER_SEC);

#endif
	bool bo = true;
	double diff = 0;
#pragma omp parallel for num_threads(8)
	for (int i = 0; i < 360; i++)
	{
		for (int j = 0; j < 180; j++)
		{
			for (int k = 0; k < 38; k++)
			{
				for (int m = 0; m < 19; m++)
				{
					bo = ((a[i][j][k][m] == b[i][j][k][m]) && bo);
					if (abs(a[i][j][k][m] - b[i][j][k][m])> diff)
					{
						printf("a[%d][%d][%d][%d]  error max\n", i, j, k, m);
						printf("a[%d][%d][%d][%d] = %g\n", i, j, k, m, a[i][j][k][m]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m, b[i][j][k][m]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m - 1, b[i][j][k][m - 1]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m + 1, b[i][j][k][m + 1]);
						printf("diff = %g\n", abs(a[i][j][k][m] - b[i][j][k][m]));
					}
					diff = max(abs(a[i][j][k][m] - b[i][j][k][m]), diff);
					//fprintf(fp3, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);
				}
			}
		}
	}
	t2 = clock();

	if (bo)
	{
		printf("Convert successfully\n");
	}
	else
	{
		printf("Max diff : %10.14lf\n", diff);
	}
	printf("Time for checking : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	return 0;
}