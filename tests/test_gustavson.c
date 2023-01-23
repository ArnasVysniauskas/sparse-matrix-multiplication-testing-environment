#include <assert.h>
#include "application/algorithms/gustavson.c"

struct ConvertToGustavsonTestCase{
    int nnz;
    int *I;
    int *J;
    double *V;
};

void test_convert_to_gustavson_matrix(){
    int no_test_cases = 3;
    int dimension = 3;

    double **test_matrices[] = create_test_matrices(no_test_cases, dimension);
    double test_matrices[3][3][3] = {
        {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}},
        {{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}},
        {{1.0, 1.0, 0.0},{0.0, 0.0, 0.0},{0.0, 1.0, 1.0}}
    };

    int test_I_1[] = {0, 0, 0, 0};
    int test_J_1[] = {};
    double test_V_1[] = {};
    int test_I_2[] = {0, 1, 2, 3};
    int test_J_2[] = {0, 1, 2};
    double test_V_2[] = {1.0, 1.0, 1.0};
    int test_I_3[] = {0, 1, 3, 4};
    int test_J_3[] = {0, 0, 2, 2};
    double test_V_3[] = {1.0, 1.0, 1.0, 1.0};

    struct ConvertToGustavsonTestCase test_cases[no_test_cases];
    test_cases[0].nnz = 0;
    test_cases[0].I = test_I_1;
    test_cases[0].J = test_J_1;
    test_cases[0].V = test_V_1;

    test_cases[1].nnz = 3;
    test_cases[1].I = test_I_2;
    test_cases[1].J = test_J_2;
    test_cases[1].V = test_V_2;

    test_cases[2].nnz = 4;
    test_cases[2].I = test_I_3;
    test_cases[2].J = test_J_3;
    test_cases[2].V = test_V_3;

    for (int t = 0; t < no_test_cases; t++)
    {  
        int *IR = safe_malloc_int(dimension + 1);
        int *JR = safe_malloc_int(test_cases[t].nnz);
        double *R = safe_malloc_double(test_cases[t].nnz);

        convert_to_gustavson_matrix(
            3, test_cases[t].nnz, test_matrices[t], IR, JR, R
        );
        for (int i = 0; i < dimension + 1; i++) assert(IR[i] == test_cases[t].I[i]);
        for (int j = 0; j < test_cases[t].nnz; j++) assert(JR[j] == test_cases[t].J[j]);
        for (int v = 0; v < test_cases[t].nnz; v++) assert(R[v] == test_cases[t].V[v]);
        
        free(IR); free(JR); free(R);
    }
}

void test_gustavson(){
    int no_test_cases = 3;
    int dimension = 3;
    double matrices_input_A[3][3][3] = {
        {{0.1, 0.0, 0.1}, {0.0, 0.1, 0.0}, {0.1, 0.0, 0.1}},
        {{0.0, 1.0, 0.0}, {2.0, 0.0, 3.0}, {4.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
    };
    int nnz_A[] = {5, 4, 0};
    double matrices_input_B[3][3][3] = {
        {{0.0, 0.1, 0.0}, {0.1, 0.0, 0.1}, {0.0, 0.1, 0.0}},
        {{0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        {{9.0, 9.0, 9.0}, {9.0, 9.0, 9.0}, {9.0, 9.0, 9.0}}
    };
    int nnz_B[] = {4, 1, 9};
    double matrices_output[3][3][3] = {
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
    };
    double matrices_excpected_output[3][3][3] = {
        {{0.0, 0.01, 0.0}, {0.02, 0.0, 0.02}, {0.0, 0.01, 0.0}},
        {{2.0, 0.0, 3.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}
    };

    for (int i = 0; i < no_test_cases; i++)
    {
        gustavson(
            dimension,
            nnz_A[i],
            nnz_B[i],
            matrices_input_A[i],
            matrices_input_B[i],
            matrices_output[i]
        );
        assert(compare_matrices(dimension, matrices_output[i], matrices_excpected_output[i]));
    }
    
}

int main(){
    printf("Starting Gustavson Test.\n");
    test_convert_to_gustavson_matrix();
    printf("'convert to gustavson matrix' tests successful.\n");
    test_gustavson();
    printf("'gustavson' tests successful\n");
}