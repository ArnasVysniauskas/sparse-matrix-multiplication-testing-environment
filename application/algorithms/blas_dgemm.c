#include <string.h>
#include <time.h>
#include <gsl/gsl_cblas.h>

#include "application/common/common.h"

void blas_dgemm(
    int dimension, int nnz_1, int nnz_2,
    double **matrix_1,
    double **matrix_2,
    struct algorithm_result *result)
{
    strcpy(result->name, __FUNCTION__);

    double *matrix_A = safe_calloc_double(dimension * dimension);
    double *matrix_B = safe_calloc_double(dimension * dimension);
    double *matrix_C = safe_calloc_double(dimension * dimension);
    for (int ci = 0; ci < dimension; ci++){
        for (int ri = 0; ri < dimension; ri++){
            matrix_A[ci * dimension + ri] = matrix_1[ci][ri];
            matrix_B[ci * dimension + ri] = matrix_2[ci][ri];
            matrix_C[ci * dimension + ri] = 0.0;   
        }
    }

    clock_t t = clock();

    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        dimension, dimension, dimension,
        (double)1.0,
        matrix_A, dimension, matrix_B, dimension,
        0.0, matrix_C, dimension
    );

    (*result).clock_ticks = (double)(clock() - t);

    for (int ci = 0; ci < dimension; ci++){
        for (int ri = 0; ri < dimension; ri++){
            result->matrix[ci][ri] = matrix_C[ci * dimension + ri];
        }
    }

    free(matrix_A);
    free(matrix_B);
    free(matrix_C);
}