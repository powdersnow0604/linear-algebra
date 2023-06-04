#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
	int *row;
	int *column;
	int rank;
}pivot;

static double** mkMatrix(size_t row, size_t column);
static void delMatrix(double** matrix);
static void delPivot(pivot* pivot);
static void printMatrix(double** matrix);

static double** rref(double** mat, pivot* piv);
static void div1row(double* matrix, size_t column, int div);
static void swap(double** a, double** b);
static void rref_sol(double** matrix);

static double** matrixMul(double** a, double** b);
static double** transpose(double** matrix);
static double** augment(double** A, double** c);
static double** reverseAugment(double** matrix, double*** c);

static void leastSquareMethod(double** matrix);

int main_linearsystem(void)
{
	//변경 시작
	size_t row;
	size_t column;
	double** matrix;
	char buf[32];

	//행렬 만들기
	FILE* file = fopen("matrix.txt", "r+t");

	fgets(buf, sizeof(buf), file);
	sscanf(buf,"%zu %zu", &row, &column);

	//공간 할당
	matrix = mkMatrix(row, column);

	//수 대입
	int line = 0;
	while (fgets(buf, sizeof(buf), file) != NULL)
	{
		char* save = (char*)buf;
		char* tok;
		int r = 0;
		while ((tok = strtok_s(NULL, " ", &save)) != NULL)
		{
			sscanf(tok, "%lf", &matrix[line][r++]);
		}
		line++;
	}
	fclose(file);

	//rref_result
	rref_sol(matrix);

	//delete
	delMatrix(matrix);

	return EXIT_SUCCESS;
}

static double** mkMatrix(size_t row, size_t column)
{
	double** matrix;
	matrix = (double**)malloc(sizeof(double*) * row);
	matrix[0] = (double*)malloc(sizeof(double) * column * row);
	for (int i = 1; i < row; i++)
	{
		matrix[i] = matrix[i - 1] + column;
	}
	return matrix;
}

static void delMatrix(double** matrix)
{
	free(matrix[0]);
	free(matrix);
}

static void delPivot(pivot* pivot)
{
	free(pivot->row);
	free(pivot->column);
}

static void printMatrix(double** matrix)
{
	size_t row = _msize(matrix) / sizeof(double*);
	size_t column = _msize(matrix[0]) / sizeof(double) / row;

	puts("\n\n");
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			printf("%lf ", matrix[i][j]);
		}
		puts("");
	}
}



static double** rref(double** mat, pivot* piv)
{
	size_t row = _msize(mat) / sizeof(double*);
	size_t column = _msize(mat[0]) / sizeof(double) /row - 1;
	pivot pivot;
	pivot.rank = 0;
	pivot.row = (int*)malloc(sizeof(int) * column);
	pivot.column = (int*)malloc(sizeof(int) * column);

	//duplicate
	double** matrix = mkMatrix(row, column + 1);
	memcpy(matrix[0], mat[0], _msize(mat[0]));

	// ~gaussion elimination
	for (int i = 0; i < row; i++)
	{

		// finding pivot
		if (pivot.rank == 0)
		{
			int j;
			for (j = 0; j < column; j++)
			{
				if (matrix[i][j] != 0)
				{
					div1row(matrix[i], column, j);
					pivot.row[pivot.rank] = i;
					pivot.column[pivot.rank] = j;
					pivot.rank += 1;
					break;
				}
				else
				{
					for (int k = i+1; k < row; k++)
					{
						if (matrix[k][j] != 0)
						{
							swap(&matrix[i], &matrix[k]);
							j--;
							break;
						}
					}
				}
			}
			if (j == column)
			{
				for (int k = i; k < row; k++)
				{
					if (matrix[k][column] != 0)
					{
						pivot.rank = -1;
						if (piv)
						{
							piv->rank = pivot.rank;
						}
						delPivot(&pivot);
						delMatrix(matrix);
						return mat;
					}
				}
				break;
			}
		}
		else
		{
			int j;
			for (j = pivot.column[pivot.rank-1]+1; j < column; j++)
			{
				if (matrix[i][j] != 0)
				{
					div1row(matrix[i], column, j);
					pivot.row[pivot.rank] = i;
					pivot.column[pivot.rank] = j;
					pivot.rank += 1;
					break;
				}
				else
				{
					for (int k = i + 1; k < row; k++)
					{
						if (matrix[k][j] != 0)
						{
							swap(&matrix[i], &matrix[k]);
							j--;
							break;
						}
					}
				}
			}
			if (j == column)
			{
				for (int k = i; k < row; k++)
				{
					if (matrix[k][column] != 0)
					{
						pivot.rank = -1;						
						if (piv)
						{
							piv->rank = pivot.rank;
						}
						delPivot(&pivot);
						delMatrix(matrix);
						return mat;
					}
				}
				break;
			}
		}


		// gaussion elemination
		for (int j = i+1; j < row; j++)
		{
			double temp = matrix[j][pivot.column[pivot.rank - 1]];
			for (int k = pivot.column[pivot.rank - 1]; k <= column; k++)
			{
				matrix[j][k] -= matrix[i][k] * temp;
			}
		}

	}
	
	// jordan elimination
	for (int i = 0; i < pivot.rank; i++)
	{
		for (int j = pivot.row[i] - 1; j >= 0; j--)
		{
			double temp = matrix[j][pivot.column[i]];
			for (int k = pivot.column[i]; k <= column; k++)
			{
				matrix[j][k] -= matrix[pivot.row[i]][k] * temp;
			}
		}
	}

	if (piv)
	{
		piv->rank = pivot.rank;
		piv->row = pivot.row;
		piv->column = pivot.column;
	}
	else
	{
		delPivot(&pivot);
	}
	return matrix;
}

static void div1row(double* matrix, size_t column, int div)
{
	double temp = matrix[div];
	for (int i = div; i <= column; i++)
	{
		matrix[i] = matrix[i] / temp;
	}

}

static void swap(double** a, double** b)
{
	double* temp = *a;
	*a = *b;
	*b = temp;
}

static void rref_sol(double** matrix)
{
	pivot pivot;
	double** rref_matrix = rref(matrix,&pivot);
	size_t row = _msize(rref_matrix) / sizeof(double*);
	size_t column = _msize(rref_matrix[0]) / sizeof(double) / row;

	if (pivot.rank == column - 1)
	{
		puts("unique solution\n");
		for (int i = 0; i < row; i++)
		{
			printf("   %lf\n", rref_matrix[i][column - 1]);
		}
		delPivot(&pivot);
		delMatrix(rref_matrix);
	}
	else if (pivot.rank == -1)
	{
		leastSquareMethod(matrix);
	}
	else
	{
		puts("infinite solution\n");
		int ppos = 0;
		for (int i = 0; i < column - 1; i++)
		{
			if (i == pivot.column[ppos])
			{
				ppos++;
				continue;
			}
			printf("x%d * [ ", i + 1);
			for (int j = 0; j < column - 1; j++)
			{


				if (j >= pivot.rank)
				{
					for (int k = j; k < column - 1; k++)
					{
						if (k == i)
							printf("%lf ", (double)1);
						else
							printf("%lf ", (double)0);
					}
					break;
				}
				else
				{
					if (j < row)
						printf("%lf ", -1 * rref_matrix[j][i]);
					else
						printf("%lf ", (double)0);
				}

			}
			puts("]T");
		}
		printf("const [ ");
		for (int j = 0; j < column - 1; j++)
		{
			if (j < row)
				printf("%lf ", rref_matrix[j][column - 1]);
			else
				printf("%lf ", (double)0);
		}
		puts("]T");
		delPivot(&pivot);
		delMatrix(rref_matrix);
	}
}



static double** matrixMul(double** a, double** b)
{
	size_t lm = _msize(a) / sizeof(double*);
	size_t ln = _msize(a[0]) / sizeof(double) / lm;

	size_t rm = _msize(b) / sizeof(double*);
	size_t rn = _msize(b[0]) / sizeof(double) / rm;

	if (ln != rm) return a;
	double** matrix = mkMatrix(lm, rn);
	for (int i = 0; i < lm; i++)
	{
		for (int j = 0; j < rn; j++)
		{
			matrix[i][j] = 0;
			for (int k = 0; k < ln; k++)
			{
				matrix[i][j] += a[i][k] * b[k][j];
			}
		}
	}

	return matrix;
}

static double** transpose(double** matrix)
{
	size_t row = _msize(matrix) / sizeof(double*);
	size_t column = _msize(matrix[0]) / sizeof(double) / row;

	double** trans = mkMatrix(column, row);

	for (int i = 0; i < column; i++)
	{
		for (int j = 0; j < row; j++)
		{
			trans[i][j] = matrix[j][i];
		}
	}

	return trans;
}

static double** augment(double** A, double** c)
{
	size_t row = _msize(A) / sizeof(double*);
	size_t column = _msize(A[0]) / sizeof(double) / row;

	double** aug = mkMatrix(row, column + 1);

	for (int i = 0; i < row; i++)
	{
		memcpy(aug[i], A[i], column * sizeof(double));
		aug[i][column] = c[i][0];
	}

	return aug;
}

static double** reverseAugment(double** matrix, double*** c)
{
	size_t row = _msize(matrix) / sizeof(double*);
	size_t column = _msize(matrix[0]) / sizeof(double) / row - 1;

	double** A = mkMatrix(row, column);
	*c = mkMatrix(row, 1);

	for (int i = 0; i < row; i++)
	{
		memcpy(A[i], matrix[i], column * sizeof(double));
		(*c)[i][0] = matrix[i][column];
	}

	return A;

}

static void leastSquareMethod(double** matrix)
{
	double** c;
	double** A = reverseAugment(matrix, &c);
	double** AT = transpose(A);

	double** ATA = matrixMul(AT, A);
	double** ATc = matrixMul(AT, c);

	delMatrix(c);
	delMatrix(A);
	delMatrix(AT);

	double** aug = augment(ATA, ATc);

	delMatrix(ATA);
	delMatrix(ATc);

	puts("[Least Square Method]");
	rref_sol(aug);

	delMatrix(aug);

}
