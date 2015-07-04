#include <stdio.h>
#include <time.h>
#include "serial_IO_para.h"

#define max(a,b) ((a>b)?a:b)
#define abs(a) (((a)>0)?(a):-(a))
int main()
{
	static double a[dim1][dim2][dim3][dim4];
	static double b[dim1][dim2][dim3][dim4];
	//fp for ascii read, fp2 for bi write, fp3 for ascii write
	FILE *fp, *fp2, *fp3;
	clock_t t1 = clock();
	if ((fp=fopen(path(A.txt),"r"))==NULL)
	{
		puts("Failed open A.txt");
		return 0;
	}
	for (int i = 0; i < dim1 * dim2 * dim3 * dim4; i++)
	{
		short int w, x, y, z;
		fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
		if (w >= dim1 || w < 0 || x >= dim2 || x < 0 || y >= dim3 || y < 0 || z >= dim4 || z < 0)
		{
			puts("Error bound\n");
			return 0;
		}
		else
		{
			fscanf(fp, "%lf\n", &(a[w][x][y][z]));
		}
		if (i%(dim2 * dim3 * dim4)==0)
		{
			printf("\r%5.2f%% Finished read", (double) i / (dim1 * dim2 * dim3 * dim4) * 100);
		}
	}
	printf("\rFinished total reading.\n");
	clock_t t2 = clock();
	printf("Time for reading A : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	fclose(fp);
	if ((fp2 = fopen(path(A_bi),"wb"))==NULL)
	{
		puts("Failed open A_bi.txt");
		return 0;
	}
	fwrite(&a, sizeof(a), 1, fp2);
	printf("Successfully convert to binary file\n");
	t1 = clock();
	printf("Time for write : %f s", (double)(t1 - t2) / CLOCKS_PER_SEC);
	fclose(fp2);
	t1 = clock();
	if ((fp3 = fopen(path(A_bi),"r"))==NULL)
	{
		puts("Failed open A_bi.txt");
		return 0;
	}
	fread(&b, sizeof(b), 1, fp3);
	fclose(fp3);
	t2 = clock();
	printf("Time for read A_bi : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	fp3 = fopen(path(A_ascii.txt), "w");
	for (int i = 0; i < dim1; i++)
	for (int j = 0; j < dim2; j++)
	for (int k = 0; k < dim3; k++)
	for (int m = 0; m < dim4; m++)
		fprintf(fp3, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);
	fclose(fp3);
	t1 = clock();
	printf("Time for write A_ascii : %f s\n", (double)(t1 - t2) / CLOCKS_PER_SEC);

	bool bo = true;
	double diff = 0;
	for (int i = 0; i < dim1; i++)
	{
		for (int j = 0; j < dim2; j++)
		{
			for (int k = 0; k < dim3; k++)
			{
				for (int m = 0; m < dim4; m++)
				{
					bo = ((a[i][j][k][m] == b[i][j][k][m]) && bo);
					if (abs(a[i][j][k][m] - b[i][j][k][m])> diff)
					{
						printf("a[%d][%d][%d][%d]  error max\n", i, j, k, m);
						printf("a[%d][%d][%d][%d] = %g\n", i, j, k, m, a[i][j][k][m]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m, b[i][j][k][m]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m-1, b[i][j][k][m-1]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m+1, b[i][j][k][m+1]);
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