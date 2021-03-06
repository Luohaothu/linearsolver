#include "parameters.h"
#include <stdio.h>
//return x value according to index
double pos_x_value(double x[dim1][dim2][dim3], int i, int j, int k, int p)
{
	switch (p + 1)
	{
	//case 1 : x(i, j, k)
	case 1: return x[i][j][k]; break;
	//case 2 : x(i - 1, j, k)
	case 2: return x[(i - 1 + dim1) % dim1][j][k]; break;
	//case 3 : x(i + 1, j, k)
	case 3: return x[(i + 1) % dim1][j][k]; break;
	//case 4 : x(i, j - 1, k)
	case 4: return j == 0 ? x[(i + dim1 / 2) % dim1][j][k] : x[i][j - 1][k]; break;
	//case 5 : x(i, j + 1, k)
	case 5: return j == dim2 - 1 ? x[(i + dim1 / 2) % dim1][j][k] : x[i][j + 1][k]; break;
	//case 6 : x(i + 1, j + 1, k)
	case 6: return j == dim2 - 1 ? x[(i + dim1 / 2 + 1) % dim1][j][k] : x[(i + 1) % dim1][j + 1][k]; break;
	//case 7 : x(i + 1, j - 1, k)
	case 7: return j == 0 ? x[(i + dim1 / 2 + 1) % dim1][j][k] : x[(i + 1) % dim1][j - 1][k]; break;
	//case 8 : x(i - 1, j - 1, k)
	case 8: return j == 0 ? x[(i + dim1 / 2 - 1) % dim1][j][k] : x[(i - 1 + dim1) % dim1][j - 1][k]; break;
	//case 9 : x(i - 1, j + 1, k)
	case 9: return j == dim2 - 1 ? x[(i + dim1 / 2 - 1) % dim1][j][k] : x[(i - 1 + dim1) % dim1][j + 1][k]; break;
	//case 10: x(i, j, k - 1)
	case 10:return k == 0 ? 0 : x[i][j][k - 1]; break;
	//case 11: x(i - 1, j, k - 1)
	case 11:return k == 0 ? 0 : x[(i - 1 + dim1) % dim1][j][k - 1]; break;
	//case 12: x(i + 1, j, k - 1)
	case 12:return k == 0 ? 0 : x[(i + 1) % dim1][j][k - 1]; break;
	//case 13: x(i, j - 1, k - 1)
	case 13:return k == 0 ? 0 : (j == 0 ? x[(i + dim1 / 2) % dim1][j][k - 1] : x[i][j - 1][k - 1]); break;
	//case 14: x(i, j + 1, k - 1)
	case 14:return k == 0 ? 0 : (j == dim2 - 1 ? x[(i + dim1 / 2) % dim1][j][k - 1] : x[i][j + 1][k - 1]); break;
	//case 15: x(i, j, k + 1)
	case 15:return k == dim3 - 1 ? 0 : x[i][j][k + 1]; break;
	//case 16: x(i - 1, j, k + 1)
	case 16:return k == dim3 - 1 ? 0 : x[(i - 1 + dim1) % dim1][j][k + 1]; break;
	//case 17: x(i + 1, j, k + 1)
	case 17:return k == dim3 - 1 ? 0 : x[(i + 1) % dim1][j][k + 1]; break;
	//case 18: x(i, j - 1, k + 1)
	case 18:return k == dim3 - 1 ? 0 : (j == 0 ? x[(i + dim1 / 2) % dim1][j][k + 1] : x[i][j - 1][k + 1]); break;
	//case 19: x(i, j + 1, k + 1)
	case 19:return k == dim3 - 1 ? 0 : (j == dim2 - 1 ? x[(i + dim1 / 2) % dim1][j][k + 1] : x[i][j + 1][k + 1]); break;
	default:return 0;
		break;
	}
}

double xvalue(double(&x)[block_step1 + 2][block_step2 + 2][dim3], int i, int j, int k, int p)
{
		int ii = i + 1;
		int jj = j + 1;
		switch (p+1)
		{
			//case 1 : x(i, j, k)
		case 1: return x[ii][jj][k]; break;
			//case 2 : x(i - 1, j, k)
		case 2: return x[ii - 1][jj][k]; break;
			//case 3 : x(i + 1, j, k)
		case 3: return x[ii + 1][jj][k]; break;
			//case 4 : x(i, j - 1, k)
		case 4: return x[ii][jj - 1][k]; break;
			//case 5 : x(i, j + 1, k)
		case 5: return x[ii][jj + 1][k]; break;
			//case 6 : x(i + 1, j + 1, k)
		case 6: return x[ii + 1][jj + 1][k]; break;
			//case 7 : x(i + 1, j - 1, k)
		case 7: return x[ii + 1][jj - 1][k]; break;
			//case 8 : x(i - 1, j - 1, k)
		case 8: return x[ii - 1][jj - 1][k]; break;
			//case 9 : x(i - 1, j + 1, k)
		case 9: return x[ii - 1][jj + 1][k]; break;
			//case 10: x(i, j, k - 1)
		case 10:return k == 0 ? 0 : x[ii][jj][k - 1]; break;
			//case 11: x(i - 1, j, k - 1)
		case 11:return k == 0 ? 0 : x[ii - 1][jj][k - 1]; break;
			//case 12: x(i + 1, j, k - 1)
		case 12:return k == 0 ? 0 : x[ii + 1][jj][k - 1]; break;
			//case 13: x(i, j - 1, k - 1)
		case 13:return k == 0 ? 0 : x[ii][jj - 1][k - 1]; break;
			//case 14: x(i, j + 1, k - 1)
		case 14:return k == 0 ? 0 : x[ii][jj + 1][k - 1]; break;
			//case 15: x(i, j, k + 1)
		case 15:return k == dim3 - 1 ? 0 : x[ii][jj][k + 1]; break;
			//case 16: x(i - 1, j, k + 1)
		case 16:return k == dim3 - 1 ? 0 : x[ii - 1][jj][k + 1]; break;
			//case 17: x(i + 1, j, k + 1)
		case 17:return k == dim3 - 1 ? 0 : x[ii + 1][jj][k + 1]; break;
			//case 18: x(i, j - 1, k + 1)
		case 18:return k == dim3 - 1 ? 0 : x[ii][jj - 1][k + 1]; break;
			//case 19: x(i, j + 1, k + 1)
		case 19:return k == dim3 - 1 ? 0 : x[ii][jj + 1][k + 1]; break;

		default:return 0;
			break;
		}
}