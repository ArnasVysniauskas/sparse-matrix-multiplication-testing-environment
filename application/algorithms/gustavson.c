#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include "application/common/common.h"

void convert_to_gustavson_matrix(
    int dimension,
    int nnz,
    double **matrix,
    int *I,
    int *J,
    double *V)
{
    int counter = 0;
    for(int ri = 0; ri < dimension; ri++)
    {
        I[ri] = counter;
        for (int ci = 0; ci < dimension; ci++)
        {
            if(compare_double(matrix[ci][ri], 0.0) == false)
            {
                J[counter] = ci;
                V[counter] = matrix[ci][ri];
                counter++;
            }
        }
    }
    I[dimension] = counter;
}

void print_gustavson_matrix(int dimension, int nnz, int *I, int *J, double *V){
    for(int i = 0; i < dimension; i ++){
        printf("%.5i ", I[i]);
        for (int j = 0; j < I[i + 1] - I[i] - 1; j++){
            printf("      ");
        }
    }
    printf("%.5i\n", I[dimension]);
    for(int i = 0; i < nnz; i ++){
        printf("%.5i ", J[i]);
    }
    printf("\n");
    for(int i = 0; i < nnz; i ++){
        printf("%5.2f ", V[i]);
    }
    printf("\n");
}

//struct algorithm_result gustavson(
void gustavson(
    int dimension, int nnz_1, int nnz_2,
    double **matrix_1,
    double **matrix_2,
    struct algorithm_result *result)
{
    strcpy((*result).name, __FUNCTION__);

    // intitial nnz_3 approximation
    int nnz_3 = min(nnz_1, dimension) * min(nnz_2, dimension);

    int *IA = safe_calloc_int(dimension + 1);
    int *JA = safe_calloc_int(nnz_1);
    double *A = safe_calloc_double(nnz_1);
    int *IB = safe_calloc_int(dimension + 1);
    int *JB = safe_calloc_int(nnz_2);
    double *B = safe_calloc_double(nnz_2);
    int *IC = safe_calloc_int(dimension + 1);
    int *JC = safe_calloc_int(nnz_3);
    double *C = safe_calloc_double(nnz_3);

    convert_to_gustavson_matrix(dimension, nnz_1, matrix_1, IA, JA, A);
    convert_to_gustavson_matrix(dimension, nnz_2, matrix_2, IB, JB, B);
    
    // Starting the algorithm
    clock_t t = clock();

    int i, ip, j, jp, k, kp, v, vp;
    int *xb = safe_calloc_int(dimension);
    double *x = safe_calloc_double(dimension);

    ip = 0;
    for (v = 0; v < dimension; v++) xb[v] = (int)-1;

    for(i = 0; i < dimension; i ++){
        IC[i] = ip;
        for(jp = IA[i]; jp < IA[i + 1]; jp++){
            j = JA[jp];
            for(kp = IB[j]; kp < IB[j + 1]; kp++){
                k = JB[kp];
                if (xb[k] != i){
                    JC[ip] = k;
                    ip = ip + 1;
                    xb[k] = i;
                    x[k] = A[jp] * B[kp];
                }else{
                    x[k] = x[k] + A[jp] * B[kp];
                }
            }
        }
        for(vp = IC[i]; vp < ip; vp++){
            v = JC[vp];
            C[vp] = x[v];
        }
    }
    IC[dimension] = ip;
    (*result).clock_ticks = (int)(clock() - t);
    
    // Building result_matrix from gustavson format
    for (int ri = 0; ri < dimension; ri++)
    {
        for (int cpi = IC[ri]; cpi < IC[ri + 1]; cpi++)
        {
            int ci = JC[cpi];
            (*result).matrix[ci][ri] = C[cpi];
        }
    }

    // freeing the memory
    free(x); x = NULL;
    free(xb); xb = NULL;
    free(C); C = NULL;
    free(JC); JC = NULL;
    free(IC); IC = NULL;
    free(B); B = NULL;
    free(JB); JB = NULL;
    free(IB); IB = NULL;
    free(A); A = NULL;
    free(JA); JA = NULL;
    free(IA); IA = NULL;
}