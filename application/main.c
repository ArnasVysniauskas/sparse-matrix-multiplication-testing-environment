#include "algorithms/gustavson.c"
#include "algorithms/naive.c"
#include "algorithms/blas_dgemm.c"
#include "algorithms/hypersparse_gemm_1.c"
#include "application/common/common.h"

bool testing = true;
int testing_quadrant_size = 100;

double sparsity_start = 0.0;
double sparsity_end = 1.0;
double sparsity_step = 0.2;

int dimension_start = 1;
int dimension_end = 10;
int dimension_step = 1;

bool activate_naive_algorithm = true;
bool activate_gustavson_algorithm = true;
bool activate_blas_dgemm_algorithm = true;
bool activate_hypersparse_gemm_algorithm_1 = true;

char results_file_name[40];
int result_id = 1;
int repeats = 1;

void environment(int dimension, double sparsity){
    //---------------------------------------------------------------------------------
    // initialising algorithms
    //---------------------------------------------------------------------------------
    int no_algorithms = 4;
    struct algorithm_result algorithms[no_algorithms];
    algorithms[0].active = activate_naive_algorithm;
    algorithms[0].matrix = create_matrix(dimension);
    algorithms[0].run = naive;
    algorithms[1].active = activate_gustavson_algorithm;
    algorithms[1].matrix = create_matrix(dimension);
    algorithms[1].run = gustavson;
    algorithms[2].active = activate_blas_dgemm_algorithm;
    algorithms[2].matrix = create_matrix(dimension);
    algorithms[2].run = blas_dgemm;
    algorithms[3].active = activate_hypersparse_gemm_algorithm_1;
    algorithms[3].matrix = create_matrix(dimension);
    algorithms[3].run = hypersparse_gemm_1;

    //---------------------------------------------------------------------------------
    // generating input matrices
    //---------------------------------------------------------------------------------
    double **matrix_1 = create_matrix(dimension);
    double **matrix_2 = create_matrix(dimension);
    int nnz_1 = matrix_generator(dimension, sparsity, matrix_1);
    int nnz_2 = matrix_generator(dimension, sparsity, matrix_2);

    //---------------------------------------------------------------------------------
    // running algorithms
    //---------------------------------------------------------------------------------
    for (int i = 0; i < no_algorithms; i++)
    {   
        if(algorithms[i].active == false) continue;
        algorithms[i].run(dimension, nnz_1, nnz_2, matrix_1, matrix_2, &(algorithms[i]));
    }

    //---------------------------------------------------------------------------------
    // Testing the first quadrant (100x100) if the results are correct
    //---------------------------------------------------------------------------------
    if(algorithms[0].active == true && testing == true)
    {
        for (int i = 1; i < no_algorithms; i++)
        {
            if(algorithms[i].active == false) continue;
            if(
                compare_matrices(
                    dimension < testing_quadrant_size ? dimension : testing_quadrant_size,
                    algorithms[0].matrix,
                    algorithms[i].matrix
                ) == true
            ) continue;

            printf(
                "Matrix result from algorithm %s doesn't match expected result\n",
                algorithms[i].name
            );
            printf("Matrix 1:\n");
            print_matrix(dimension, matrix_1);
            printf("Matrix 2:\n");
            print_matrix(dimension, matrix_2);
            printf("Expected result:\n");
            print_matrix(dimension, algorithms[0].matrix);
            printf("Result:\n");
            print_matrix(dimension, algorithms[i].matrix);
            
        }
    }
    //---------------------------------------------------------------------------------
    // Saving results
    //---------------------------------------------------------------------------------
    
    FILE *f = fopen(results_file_name, "a");
    if (f == NULL)
    {
        printf("error opening file\n");
    } else 
    {
        fprintf(f, "%d, %f", dimension, sparsity);
        for (int i = 0; i < no_algorithms; i++)
        {
            if(algorithms[i].active == false) continue;
            fprintf(f, ", %d", algorithms[i].clock_ticks);
        }
        fprintf(f, "\n");
        fclose(f);
    }

    //---------------------------------------------------------------------------------
    // Deallocating memory
    //---------------------------------------------------------------------------------
    for (int i = 0; i < no_algorithms; i++) delete_matrix(dimension, algorithms[i].matrix);
    delete_matrix(dimension, matrix_1);
    delete_matrix(dimension, matrix_2);
}

double get_config_value(char *key, FILE *config){
    char buffer[50];
    char value_s[50];
    if(fgets(buffer, 50, config) == NULL){
        key[0] = '\0';
        return 0.0;
    }

    int key_index = 0;
    char key_char = buffer[key_index];
    while(key_char != '='){
        key[key_index] = key_char;
        key_index++;
        key_char = buffer[key_index];
    }
    key[key_index] = '\0';
    key_index++;
    int key_length = key_index;
    key_char = buffer[key_index];
    while(key_char != '\n'){
        value_s[key_index - key_length] = key_char;
        key_index++;
        key_char = buffer[key_index];
    }
    value_s[key_index - key_length] = '\0';
    return strtof(value_s, NULL);
}

void read_config(){
    FILE *config = fopen("config.txt", "r");
    if (config == NULL)
    {
        printf("error opening config file\n");
    } else 
    {
        char *key_ptr = calloc(50, sizeof(char));
        double value = get_config_value(key_ptr, config);

        while(key_ptr[0] != '\0'){
            if(strcmp(key_ptr, "testing") == 0) testing = (bool)value;
            else if(strcmp(key_ptr, "testing_quadrant_size") == 0) testing_quadrant_size = (int)value;
            else if(strcmp(key_ptr, "sparsity_start") == 0) sparsity_start = (double)value;
            else if(strcmp(key_ptr, "sparsity_end") == 0) sparsity_end = (double)value;
            else if(strcmp(key_ptr, "sparsity_step") == 0) sparsity_step = (double)value;
            else if(strcmp(key_ptr, "dimension_start") == 0) dimension_start = (int)value;
            else if(strcmp(key_ptr, "dimension_end") == 0) dimension_end = (int)value;
            else if(strcmp(key_ptr, "dimension_step") == 0) dimension_step = (int)value;
            else if(strcmp(key_ptr, "naive_algorithm") == 0) activate_naive_algorithm = (bool)value;
            else if(strcmp(key_ptr, "gustavson_algorithm") == 0) activate_gustavson_algorithm = (bool)value;
            else if(strcmp(key_ptr, "dgemm_algorithm") == 0) activate_blas_dgemm_algorithm = (bool)value;
            else if(strcmp(key_ptr, "hypersparse_gemm_algorithm_1") == 0) activate_hypersparse_gemm_algorithm_1 = (bool)value;
            else if(strcmp(key_ptr, "result_id") == 0) result_id = (int)value;
            else if(strcmp(key_ptr, "repeats") == 0) repeats = (int)value;
            else printf("Unexpected key in config. Key: '%s'\n", key_ptr);

            value = get_config_value(key_ptr, config);
        }

        fclose(config);
        free(key_ptr);
    }
}

int main(int argc, const char *argv[]){
    read_config();
    sprintf(results_file_name, "results_%d.csv", result_id);

    FILE *f = fopen(results_file_name, "a");
    if (f == NULL)
    {
        printf("error opening file\n");
    } else 
    {
        fprintf(f, "dimension, sparsity");
        if (activate_naive_algorithm == true) fprintf(f, ", naive");
        if (activate_gustavson_algorithm == true) fprintf(f, ", gustavson");
        if (activate_blas_dgemm_algorithm == true) fprintf(f, ", blas_dgemm");
        if (activate_hypersparse_gemm_algorithm_1 == true) fprintf(f, ", hypersparse_gemm_1");
        fprintf(f, "\n");
        fclose(f);
    }
    for (int r = 0; r < repeats; r++){
        for (double sparsity = sparsity_start; sparsity < sparsity_end; sparsity += sparsity_step){
            for (int i = dimension_start; i < dimension_end; i += dimension_step){
                printf("%d) Test case -> dimension = %d, sparsity = %.5f\n", r, i, sparsity);
                environment(i, sparsity);
            }
        }
    }

    return 0;
}