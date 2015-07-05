#include "IO_test.h"
int convert()
{
	static double a[dim1][dim2][dim3][dim4];
	static double b[dim1][dim2][dim3];
	static double x0[dim1][dim2][dim3];
	//file pointer for read ascii data
	FILE *r_fp;
	//file pointer for write binary data
	FILE *w_fp1, *w_fp2, *w_fp3, *w_fp4;
	//read A
	if ((r_fp = fopen(path(A.txt), "r")) == NULL)
	{
		puts("Failed open A.txt");
		return 10;
	}
	for (int i = 0; i < dim1 * dim2 * dim3 * dim4; i++)
	{
		short int w, x, y, z;
		fscanf(r_fp, "%d%d%d%d", &w, &x, &y, &z);
		if (w >= dim1 || w < 0 || x >= dim2 || x < 0 || y >= dim3 || y < 0 || z >= dim4 || z < 0)
		{
			puts("Error bound\n");
			return 11;
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
	fclose(r_fp);
	//read b
	if ((r_fp = fopen(path(b.txt), "r")) == NULL)
	{
		puts("Failed open b.txt");
		return 20;
	}
	for (int i = 0; i < dim1 * dim2 * dim3; i++)
	{
		short int x, y, z;
		fscanf(r_fp, "%d%d%d", &x, &y, &z);
		if (x >= dim1 || x < 0 || y >= dim2 || y < 0 || z >= dim3 || z < 0)
		{
			puts("Error bound\n");
			return 21;
		}
		else
		{
			fscanf(r_fp, "%lf\n", &(b[x][y][z]));
		}
		if (i % (dim2 * dim3) == 0)
		{
			printf("\r%5.2f%% Finished read", (double)i / (dim1 * dim2 * dim3) * 100);
		}
	}
	printf("\rFinished total reading.\n");
	fclose(r_fp);
	//read x0
	if ((r_fp = fopen(path(x0.txt), "r")) == NULL)
	{
		puts("Failed open x0.txt");
		return 30;
	}
	for (int i = 0; i < dim1 * dim2 * dim3; i++)
	{
		short int x, y, z;
		fscanf(r_fp, "%d%d%d", &x, &y, &z);
		if (x >= dim1 || x < 0 || y >= dim2 || y < 0 || z >= dim3 || z < 0)
		{
			puts("Error bound\n");
			return 31;
		}
		else
		{
			fscanf(r_fp, "%lf\n", &(x0[x][y][z]));
		}
		if (i % (dim2 * dim3) == 0)
		{
			printf("\r%5.2f%% Finished read", (double)i / (dim1 * dim2 * dim3) * 100);
		}
	}
	printf("\rFinished total reading.\n");
	fclose(r_fp);
	//write A_bi
	if (   (w_fp1 = fopen(path(A_bi.part1), "wb")) == NULL
		|| (w_fp2 = fopen(path(A_bi.part2), "wb")) == NULL
		|| (w_fp3 = fopen(path(A_bi.part3), "wb")) == NULL
		|| (w_fp4 = fopen(path(A_bi.part4), "wb")) == NULL)
	{
		puts("Failed open A_bi");
		return -1;
	}
	fwrite(&a, sizeof(a) / 4, 1, w_fp1);
	fwrite(&(a[dim1 / 4]), sizeof(a) / 4, 1, w_fp2);
	fwrite(&(a[dim1 / 2]), sizeof(a) / 4, 1, w_fp3);
	fwrite(&(a[3 * dim1 / 4]), sizeof(a) / 4, 1, w_fp4);
	printf("Successfully convert to binary file\n");
	fclose(w_fp1);
	fclose(w_fp2);
	fclose(w_fp3);
	fclose(w_fp4);
	//write b_bi
	if (   (w_fp1 = fopen(path(b_bi.part1), "wb")) == NULL
		|| (w_fp2 = fopen(path(b_bi.part2), "wb")) == NULL
		|| (w_fp3 = fopen(path(b_bi.part3), "wb")) == NULL
		|| (w_fp4 = fopen(path(b_bi.part4), "wb")) == NULL)
	{
		puts("Failed open b_bi");
		return -2;
	}
	fwrite(&b, sizeof(b) / 4, 1, w_fp1);
	fwrite(&(b[dim1 / 4]), sizeof(b) / 4, 1, w_fp2);
	fwrite(&(b[dim1 / 2]), sizeof(b) / 4, 1, w_fp3);
	fwrite(&(b[3 * dim1 / 4]), sizeof(b) / 4, 1, w_fp4);
	printf("Successfully convert to binary file\n");
	fclose(w_fp1);
	fclose(w_fp2);
	fclose(w_fp3);
	fclose(w_fp4);
	//write x0_bi
	if (   (w_fp1 = fopen(path(x0_bi.part1), "wb")) == NULL
		|| (w_fp2 = fopen(path(x0_bi.part2), "wb")) == NULL
		|| (w_fp3 = fopen(path(x0_bi.part3), "wb")) == NULL
		|| (w_fp4 = fopen(path(x0_bi.part4), "wb")) == NULL)
	{
		puts("Failed open x0_bi");
		return -3;
	}
	fwrite(&x0, sizeof(x0) / 4, 1, w_fp1);
	fwrite(&(x0[dim1 / 4]), sizeof(x0) / 4, 1, w_fp2);
	fwrite(&(x0[dim1 / 2]), sizeof(x0) / 4, 1, w_fp3);
	fwrite(&(x0[3 * dim1 / 4]), sizeof(x0) / 4, 1, w_fp4);
	printf("Successfully convert to binary file\n");
	fclose(w_fp1);
	fclose(w_fp2);
	fclose(w_fp3);
	fclose(w_fp4);

	return 0;
}