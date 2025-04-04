#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** allocate_matrix(int size) 
{
    double** matrix = (double**)malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) 
    {
        matrix[i] = (double*)malloc(size * sizeof(double));
    }
    return matrix;
}

void free_matrix(double** matrix, int size) 
{
    for (int i = 0; i < size; i++) 
    {
        free(matrix[i]);
    }
    free(matrix);
}

void mat_mult(int R1, int C1, double ** m1, int R2, int C2, double m2[][C2], double result[][C2])
{
	for (int i = 0; i < R1; i++) 
	{
		for (int j = 0; j < C2; j++) 
		{
			result[i][j] = 0;
			for (int k = 0; k < R2; k++) 
			{
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

void mat_mult1(int R1, int C1, double m1[][4], int R2, int C2, double m2[][C2], double result[][C2])
{
	for (int i = 0; i < R1; i++) 
	{
		for (int j = 0; j < C2; j++) 
		{
			result[i][j] = 0;
			for (int k = 0; k < R2; k++) 
			{
				result[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}

double determinant(int size, double** matrix) 
{
    if (size == 1) return matrix[0][0];

    double det = 0.0;
    double** minor = allocate_matrix(size - 1);
    
    for (int col = 0; col < size; col++) 
    {
        int m = 0, n = 0;
        for (int i = 1; i < size; i++) 
        {
            for (int j = 0; j < size; j++) 
            {
                if (j == col) continue;
                minor[m][n++] = matrix[i][j];
                if (n == size - 1) 
                {
                    n = 0;
                    m++;
                }
            }
        }
        det += (col % 2 == 0 ? 1 : -1) * matrix[0][col] * determinant(size - 1, minor);
    }
    
    free_matrix(minor, size - 1);
    return det;
}

void cofactor(int size, double** matrix, double** cofactors) 
{
    double** minor = allocate_matrix(size - 1);
    
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            int m = 0, n = 0;
            for (int x = 0; x < size; x++) 
            {
                for (int y = 0; y < size; y++) 
                {
                    if (x != i && y != j) 
                    {
                        minor[m][n++] = matrix[x][y];
                        if (n == size - 1) 
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            cofactors[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * determinant(size - 1, minor);
        }
    }
    
    free_matrix(minor, size - 1);
}

void transpose(int size, double** matrix, double** result) {
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            result[j][i] = matrix[i][j];
        }
    }
}

void inverse_matrix(int size, double** matrix) 
{
    double det = determinant(size, matrix);
    if (fabs(det) < 1e-9) 
    {
        printf("No Inverse\n");
        return;
    }

    double** cofactors = allocate_matrix(size);
    double** adjugate = allocate_matrix(size);
    double** inverse = allocate_matrix(size);

    cofactor(size, matrix, cofactors);
    transpose(size, cofactors, adjugate);
    
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            inverse[i][j] = adjugate[i][j] / det;
        }
    }
    
    free_matrix(cofactors, size);
    free_matrix(adjugate, size);
    free_matrix(inverse, size);
}

double distance(double x1, double y1, double x2, double y2) 
{
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

int main() 
{
    int n;
    scanf("%d", &n);
    double P[n][2];
    for(int i=0; i<n; ++i)
    {
        scanf("%lf%lf", &P[i][0], &P[i][1]);
    }
    double P1[2][2];
    for(int i=0; i<2; ++i)
    {
        scanf("%lf%lf", &P1[i][0], &P1[i][1]);
    }
    double t[n];
    t[0] = 0;
    for(int i=1; i<n; ++i)
    {
        t[i] = distance(P[i][0], P[i][1], P[i-1][0], P[i-1][1]);
    }
    double T[n][2];
    double** M = allocate_matrix(n);
    double R[n][2];
    for(int i=0; i<n; ++i)
    {
        M[0][i] = M[n-1][i] = 0;
    }
    M[0][0] = M[n-1][n-1] = 1;
    for(int i=1; i<n-1; ++i)
    {
        for(int j=0; j<n; ++j)
        {
            M[i][j] = 0;
            if(i==j+1)
            {
                M[i][j] = t[i+1];
            }
            else if(i==j)
            {
                M[i][j] = 2*(t[i+1] + t[i]);
            }
            else if(i==j-1)
            {
                M[i][j] = t[i];
            }
        }
    }
    R[0][0] = P1[0][0];
    R[0][1] = P1[0][1];
    R[n-1][0] = P1[1][0];
    R[n-1][1] = P1[1][1];

    for(int i=1; i<n-1; ++i)
    {
        R[i][0] = 3 * ((t[i]*t[i]*(P[i+1][0]-P[i][0])) + (t[i+1]*t[i+1]*(P[i][0]-P[i-1][0]))) / (t[i] * t[i+1]);
        R[i][1] = 3 * ((t[i]*t[i]*(P[i+1][1]-P[i][1])) + (t[i+1]*t[i+1]*(P[i][1]-P[i-1][1]))) / (t[i] * t[i+1]);
    }

    inverse_matrix(n, M);
    mat_mult(n, n, M, n, 2, R, T);
    
    for(int i=0; i<n-1; ++i)
    {
        printf("Segment %d: \n", i+1);
        double lol = 1 / t[i+1];
        double B1[4][4] = {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {-3 * pow(lol, 2), -2 * lol, 3 * pow(lol, 2), -1 * lol},
            {2 * pow(lol, 3), pow(lol, 2), -2 * pow(lol, 3), pow(lol, 2)}
        };
        double B2[4][2] = {
            {P[i][0], P[i][1]},
            {T[i][0], T[i][1]},
            {P[i+1][0], P[i+1][1]},
            {T[i+1][0], T[i+1][1]}
        };
        double B[4][2];
        mat_mult1(4, 4, B1, 4, 2, B2, B);
        for(double t1=0; t1<=t[i+1]; t1+=0.33)
        {
            double again[1][4] = {1, t1, t1*t1, t1*t1*t1};
            double my_pt[1][2];
            mat_mult1(1, 4, again, 4, 4, B, my_pt);
            printf("\tt=%.2lf\t%.2lf %.2lf\n", t1, my_pt[0][0], my_pt[0][1]);
        }
        printf("\n");
    }
    free_matrix(M, n);
    return 0;
}

