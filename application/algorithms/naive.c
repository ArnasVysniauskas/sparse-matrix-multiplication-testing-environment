#include <string.h>
#include <time.h>
#include "application/common/common.h"

//struct algorithm_result naive(
void naive(
    int dimension,
    int nnz_1,
    int nnz_2,
    double **matrix_1,
    double **matrix_2,
    struct algorithm_result *result
){
    strcpy((*result).name, __FUNCTION__);

    clock_t t = clock();
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            for (int k = 0; k < dimension; k++)
            {
                (*result).matrix[i][j] +=  matrix_1[k][j] * matrix_2[i][k];
            }
        }
    }
    (*result).clock_ticks = (int)(clock() - t);
    //return result;
}