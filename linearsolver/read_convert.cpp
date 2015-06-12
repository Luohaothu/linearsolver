#include <stdio.h>
#include <time.h>


#define max(a,b) ((a>b)?a:b)
#define abs(a) (((a)>0)?(a):-(a))
int main()
{
	static double a[360][180][38][19];
	static double b[360][180][38][19];
	FILE *fp, *fp2, *fp3;
	clock_t t1 = clock();
	if ((fp=fopen("D:\\linear\\linearSolverData\\case1\\A.txt","r"))==NULL)
	{
		puts("Failed open A.txt");
		return 0;
	}
	for (int i = 0; i < 360*180*38*19; i++)
	{
		short int w, x, y, z;
		fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
		if (w >= 360 || w < 0 || x >= 180 || x < 0 || y >= 38 || y < 0 || z >= 19 || z < 0)
		{
			puts("Error bound\n");
			return 0;
		}
		else
		{
			fscanf(fp, "%lf\n", &(a[w][x][y][z]));
			//printf("%lf\n", a[w][x][y][z]);
		}
		if (i%(10*180*38*19)==0)
		{
			printf("%d Finished read\n", i / (10 * 180 * 38 * 19));
		}
	}
	printf("Finished total reading.\n");
	clock_t t2 = clock();
	printf("Time for reading A : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	fclose(fp);
	if ((fp2 = fopen("D:\\linear\\linearSolverData\\case1\\A_bi.txt","wb"))==NULL)
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
	if ((fp3 = fopen("D:\\linear\\linearSolverData\\case1\\A_bi.txt","r"))==NULL)
	{
		puts("Failed open A_bi.txt");
		return 0;
	}
	fread(&b, sizeof(b), 1, fp3);
	fclose(fp3);
	t2 = clock();
	printf("Time for read A_bi : %f s\n", (double)(t2 - t1) / CLOCKS_PER_SEC);
	fp3 = fopen("D:\\linear\\linearSolverData\\case1\\A_ascii.txt", "w");

	bool bo = true;
	double diff = 0;
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
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m-1, b[i][j][k][m-1]);
						printf("b[%d][%d][%d][%d] = %g\n", i, j, k, m+1, b[i][j][k][m+1]);
						printf("diff = %g\n", abs(a[i][j][k][m] - b[i][j][k][m]));
					}
					diff = max(abs(a[i][j][k][m] - b[i][j][k][m]), diff);
					fprintf(fp3, "%hd   %hd   %hd   %hd   %2.30lf\n", i, j, k, m, b[i][j][k][m]);
				}
			}
		}
	}
	t1 = clock();
	fclose(fp3);
	
	if (bo)
	{
		printf("Convert successfully\n");
	}
	else
	{
		printf("Max diff : %10.14lf\n", diff);
	}
		printf("Time for checking : %f s\n", (double)(t1 - t2) / CLOCKS_PER_SEC);
	return 0;
}