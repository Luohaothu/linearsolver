#include "IO_test.h"

//return x value according to index
double xvalue(double x[dim1][dim2][dim3], int i, int j, int k, int p)
{
	switch (p)
	{
	case 1: return x[i][j][k]; break;
	case 2: return x[(i - 1) % dim1][j][k]; break;
	case 3: return x[(i + 1) % dim1][j][k]; break;
	case 4: return j == 0 ? x[(i - dim1 / 2) % dim1][j][k] : x[i][j - 1][k]; break;
	case 5: return j == dim1 - 1 ? x[(i - dim1 / 2) % dim1][j][k] : x[i][j + 1][k]; break;
	case 6: return j == dim1 - 1 ? 0 : x[(i + 1) % dim1][j + 1][k]; break;
	case 7: return j == 0 ? 0 : x[(i + 1) % dim1][j - 1][k]; break;
	case 8: return j == 0 ? 0 : x[(i - 1) % dim1][j - 1][k]; break;
	case 9: return j == dim1 - 1 ? 0 : x[(i - 1) % dim1][j + 1][k]; break;
	case 10:return k == 0 ? 0 : x[i][j][k - 1]; break;
	case 11:return k == 0 ? 0 : x[(i - 1) % dim1][j][k - 1]; break;
	case 12:return k == 0 ? 0 : x[(i + 1) % dim1][j][k - 1]; break;
	case 13:return k == 0 ? 0 : (j == 0 ? x[(i - dim1 / 2) % dim1][j][k - 1] : x[i][j - 1][k - 1]); break;
	case 14:return k == 0 ? 0 : (j == dim1 - 1 ? x[(i - dim1 / 2) % dim1][j][k] : x[i][j + 1][k - 1]); break;
	case 15:return k == dim1 - 1 ? 0 : x[i][j][k + 1]; break;
	case 16:return k == dim1 - 1 ? 0 : x[(i - 1) % dim1][j][k + 1]; break;
	case 17:return k == dim1 - 1 ? 0 : x[(i + 1) % dim1][j][k + 1]; break;
	case 18:return k == dim1 - 1 ? 0 : (j == 0 ? x[(i - dim1 / 2) % dim1][j][k + 1] : x[i][j - 1][k + 1]); break;
	case 19:return k == dim1 - 1 ? 0 : (j == dim1 - 1 ? x[(i - dim1 / 2) % dim1][j][k + 1] : x[i][j + 1][k + 1]); break;
	default:return 0;
		break;
	}
}