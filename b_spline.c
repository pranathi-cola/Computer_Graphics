#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N_CONTROL 4         
#define DEGREE 3 
#define NUM_POINTS 200 

double P[N_CONTROL][2] = {
    {0.0, 0.0},
    {0.0, 1.0},
    {1.0, 1.0},
    {1.0, 0.0}
};

double U[N_CONTROL + DEGREE + 1] = {0, 0, 0, 0, 1, 1, 1, 1};

double N_basis(int i, int p, double t, double* U) {
    if (p == 0) {
        if ((U[i] <= t && t < U[i+1]) || (t == U[N_CONTROL + DEGREE] && t == U[i+1]))
            return 1.0;
        else
            return 0.0;
    } else {
        double coef1 = U[i+p] - U[i];
        double coef2 = U[i+p+1] - U[i+1];

        double term1 = 0.0, term2 = 0.0;

        if (coef1 != 0)
            term1 = ((t - U[i]) / coef1) * N_basis(i, p-1, t, U);
        if (coef2 != 0)
            term2 = ((U[i+p+1] - t) / coef2) * N_basis(i+1, p-1, t, U);

        return term1 + term2;
    }
}

int main() 
{
    double t_start = U[DEGREE];
    double t_end = U[N_CONTROL];  

    for (int j = 0; j < NUM_POINTS; j++) {
        double t = t_start + (t_end - t_start) * j / (NUM_POINTS - 1);
        double C[2] = {0.0, 0.0};

        for (int i = 0; i < N_CONTROL; i++) {
            double Ni = N_basis(i, DEGREE, t, U);
            C[0] += Ni * P[i][0];
            C[1] += Ni * P[i][1];
        }

        printf("%lf,%lf\n", C[0], C[1]);  
    }

    return 0;
}
