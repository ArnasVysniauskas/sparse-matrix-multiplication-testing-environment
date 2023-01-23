#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <stdlib.h>
#include <stdio.h> 
#include <stdbool.h>
#include <time.h>

#define CLOCKS_PER_NS CLOCKS_PER_SEC / 1000000
#define CLOCKS_PER_MS CLOCKS_PER_SEC / 1000

int *safe_calloc_int(int size){
    void *p = calloc(size, sizeof(int));
    if (p != NULL) return p;
    printf(
        "Fatal: failed to allocate %zu bytes.\n",
        size * sizeof(int)
    );
    abort();
}

int *safe_realloc_int(int *p, int new_size){
    p = realloc(p, new_size * sizeof(int));
    if(p != NULL) return p;
    printf(
        "FATAL: failed to reallocate %zu bytes.\n",
        new_size * sizeof(double)
    );
    abort();
}

void *safe_calloc(int count, size_t size){
    void *p = calloc(count, size);
    if (p != NULL) return p;
    printf(
        "Fatal: failed to allocate %zu bytes.\n",
        count * size
    );
    abort();
}

double *safe_calloc_double(int size){
    void *p = calloc(size, sizeof(double));
    if (p != NULL) return p;
    printf(
        "Fatal: failed to allocate %zu bytes.\n",
        size * sizeof(double)
    );
    abort();    
}

void *safe_realloc_double(void *p, int new_size){
    p = realloc(p, new_size * sizeof(double));
    if(p != NULL) return p;
    printf(
        "FATAL: failed to reallocate %zu bytes.\n",
        new_size * sizeof(double)
    );
    abort();
}

bool compare_double(double a, double b){
    float epsilon = 1e-6;
    if(b - epsilon < a && a < b + epsilon) return true;
    return false;
}

bool compare_matrices(int dimension, double **matrix_1, double **matrix_2){
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            if(compare_double(matrix_1[i][j], matrix_2[i][j]) == false) return false;
        }
    }
    return true;
}

double rand_double(){
    int precision = 1000000;
    if(precision > RAND_MAX){
        printf("Precision of %d is invalid, max precision %d\n", precision, RAND_MAX);
        abort();
    }
    return (double)(rand() % precision + 1) / precision;
}

int rand_int(int start, int end){
    // returns from [start; end)
    if(end - start > RAND_MAX){
        printf("The difference end - start = %d is too big, maximum difference %d\n", end - start, RAND_MAX);
        abort();
    }
    return (rand() % (end - start)) + start;
}

bool rand_bool(double true_weigth){
    return true_weigth > rand_double() ? true : false;
}

int min(int a, int b){
    return a < b? a : b;
}

void delete_matrix(int dimension, double **matrix_buffer){
    for (int i = 0; i < dimension; i++) free(matrix_buffer[i]);
    free(matrix_buffer);
}

void reset_matrix(int dimension, double **matrix){
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            matrix[i][j] = 0.0;
        }
    }
}

double **create_matrix(int dimension){
    double **matrix_buffer = malloc(dimension * sizeof(double *));
    if(matrix_buffer == NULL)
    {
        printf("failed allocating memory for matrix\n");
        exit(1);    
    }
    for (int i = 0; i < dimension; i++)
    {
        matrix_buffer[i] = malloc(dimension * sizeof(double));
        if (matrix_buffer[i] == NULL)
        {
            printf("failed allocating memory for matrix\n");
            delete_matrix(i, matrix_buffer);
            exit(1);
        } 
    }
    reset_matrix(dimension, matrix_buffer);
    return matrix_buffer;
}

int matrix_generator(int dimension, double sparsity, double **matrix){
    clock_t start = clock();
    int max_index = dimension * dimension;
    int nnz = compare_double(sparsity, 0.0) == true? max_index: (int)(max_index * (1.0 - sparsity)) + 1;
    int *free_spaces = malloc((max_index - 1) * sizeof(int));
    for (int i = 0; i < max_index; i++) free_spaces[i] = i;
    
    for (int i = 0; i < nnz; i++)
    {
        int ptr_index = rand_int(0, max_index);
        int index = free_spaces[ptr_index];
        int row_i = index % dimension;
        int column_i = (int)((index - row_i) / dimension);

        matrix[row_i][column_i] = rand_double();
        free_spaces[ptr_index] = free_spaces[max_index - 1];
        max_index -= 1;
    }
    printf("time = %d\n", (int)(clock() - start));
    return nnz;
}

void print_matrix(int dimension, double **matrix){
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            printf("%.4f ", matrix[j][i]);
        }
        printf("\n");
    }
    printf("\n");
}

void copy_int_vector(int *destination, int *source, int size){
    for(int i = 0; i < size; i++){
        destination[i] = source[i];
    }
}

void copy_double_vector(double *destination, double *source, int size){
    for(int i = 0; i < size; i++){
        destination[i] = source[i];
    }
}

struct algorithm_result{
    bool active;
    char name[20];
    int clock_ticks;
    double **matrix;
    void (*run)(int, int, int, double**, double**, struct algorithm_result*);
};

#endif