
#include <string.h>
#include <stdbool.h>
#include "application/common/common.h"
#include "application/algorithms/algorithms_common.c"

struct Intersection{
    int index;
    int sizeA;
    int sizeB;
    int startA;
    int startB;
    int curA;
    int curB;
    bool empty;
};

void print_intersections(struct Intersection *intersections, int intersections_count){
    printf("i, szA, szB, stA, stB, crA, crB, empty\n");
    for(int i = 0; i < intersections_count; i++){
        printf("%.1i, %.3i, %.3i, %.3i, %.3i, %.3i, %.3i, %.3i\n",
            intersections[i].index,
            intersections[i].sizeA,
            intersections[i].sizeB,
            intersections[i].startA,
            intersections[i].startB,
            intersections[i].curA,
            intersections[i].curB,
            intersections[i].empty);
    }
    printf("\n");
}

void print_csm_matrix(int dimension, int nnz, int *I, int *J, double *V){
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

void print_dcsm_matrix(struct DCSMatrix *dcsm, int nnz){
    int aux_i = 0;
    do{
        printf("%d, ", dcsm->AUX[aux_i]);
        aux_i++;
    } while(dcsm->AUX[aux_i] != dcsm->AUX[aux_i - 1]);
    printf("%d\n", dcsm->AUX[aux_i]);
    for(int i = 0; i < aux_i - 1; i++){
        printf("%d, ", dcsm->J[i]);
    }
    printf("%d\n", dcsm->J[aux_i - 1]);
    for(int i = 0; i < aux_i; i++){
        printf("%d, ", dcsm->P[i]);
    }
    printf("%d\n", dcsm->P[aux_i]);
    for(int i = 0; i < nnz - 1; i++){
        printf("%d, ", dcsm->I[i]);
    }
    printf("%d\n", dcsm->I[nnz - 1]);
    for(int i = 0; i < nnz - 1; i++){
        printf("%.3f, ", dcsm->NUM[i]);
    }
    printf("%.3f\n", dcsm->NUM[nnz - 1]);
}

int find_intersections(
    struct Intersection *intersections,
    struct DCSMatrix *dcsc,
    struct DCSMatrix *dcsr,
    int nnz_1, int nnz_2
){
    int nzc, nzr;
    int intersections_count = 0;
    int i = 0;
    while(dcsc->AUX[i] != dcsc->AUX[i + 1]) i++;
    nzc = dcsc->AUX[i];
    i = 0;
    while(dcsr->AUX[i] != dcsr->AUX[i + 1]) i++;
    nzr = dcsr->AUX[i];
    int idx_1 = 0;
    int idx_2 = 0;
    while(idx_1 < nzc && idx_2 < nzr){
        if(dcsc->J[idx_1] == dcsr->J[idx_2]){
            intersections[intersections_count].index = dcsc->J[idx_1];
            intersections[intersections_count].sizeA = dcsc->P[idx_1 + 1] - dcsc->P[idx_1];
            intersections[intersections_count].sizeB = dcsr->P[idx_2 + 1] - dcsr->P[idx_2];
            intersections[intersections_count].startA = dcsc->P[idx_1];
            intersections[intersections_count].startB = dcsr->P[idx_2];
            intersections[intersections_count].curA = dcsc->P[idx_1];
            intersections[intersections_count].curB = dcsr->P[idx_2];
            intersections[intersections_count].empty = false;
            idx_1++;
            idx_2++;
            intersections_count++;
        }else if(dcsc->J[idx_1] < dcsr->J[idx_2]){
            idx_1++;
        }else{
            idx_2++;
        }
    }
    printf("nzc=%d, nzr=%d\n", nzc, nzr);
    return intersections_count;
}

void mult_intersection_1(
    struct DCSMatrix *dcsc_1,
    struct DCSMatrix *dcsc_2,
    struct Intersection *intersections,
    struct HeapNode *heap,
    int *heap_size,
    int *ti
){
    //printf("TI TEST 2: ti = %d\n", *ti);
    int ptr_a = (*dcsc_1).P[intersections[*ti].curA];
    int ptr_b = (*dcsc_2).P[intersections[*ti].curB];
    struct HeapNode *new_node = malloc(sizeof(struct HeapNode));
    //printf("TI TEST 2.1: ti = %d\n", *ti);
    new_node->key[0] = (*dcsc_1).I[ptr_a];
    new_node->key[1] = (*dcsc_2).J[ptr_b];
    new_node->value = (*dcsc_1).NUM[ptr_a] * (*dcsc_2).NUM[ptr_b];
    new_node->index = *ti;
    printf("TI TEST 2.2: ti = %d\n", *ti);
    insert(heap, heap_size, new_node, ti);
    free(new_node); new_node = NULL;
    //printf("TI TEST 3: ti = %d\n", *ti);
    //printf("TI TEST 3.1: ti = %d\n", *ti);
}

void increment_list(struct Intersection *intersections, int ti){
    printf("TEST: increment_list at index %d\n", ti);
    if(intersections[ti].curB - intersections[ti].startB < intersections[ti].sizeB){
        intersections[ti].curB = intersections[ti].curB + 1;
    }else{
        intersections[ti].curB = intersections[ti].startB;
        if(intersections[ti].curA - intersections[ti].startA < intersections[ti].sizeA){
            intersections[ti].curA = intersections[ti].curA + 1;
        }else{
            intersections[ti].empty = true;
        }
    }
}

bool intersections_empty(struct Intersection *intersections, int intersections_count){
    for(int i = 0; i < intersections_count; i ++){
        if(intersections[i].empty == false) return false;
    }
    return true;
}

void hypersparse_gemm_1(
    int dimension, int nnz_1, int nnz_2,
    double **matrix_1,
    double **matrix_2,
    struct algorithm_result *result)
{
    strcpy((*result).name, __FUNCTION__);
    struct CSMatrix *csc_1 = malloc(sizeof(struct CSMatrix));
    struct CSMatrix *csr_1 = malloc(sizeof(struct CSMatrix));
    struct CSMatrix *csc_2 = malloc(sizeof(struct CSMatrix));
    struct CSMatrix *csr_2 = malloc(sizeof(struct CSMatrix));
    struct DCSMatrix *dcsc_1 = malloc(sizeof(struct DCSMatrix));
    struct DCSMatrix *dcsr_1 = malloc(sizeof(struct DCSMatrix));
    struct DCSMatrix *dcsc_2 = malloc(sizeof(struct DCSMatrix));
    struct DCSMatrix *dcsr_2 = malloc(sizeof(struct DCSMatrix));

    init_csm(csc_1, dimension, nnz_1);
    convert_to_csc(dimension, nnz_1, matrix_1, csc_1);
    init_csm(csr_1, dimension, nnz_1);
    convert_to_csr(dimension, nnz_1, matrix_1, csr_1);
    init_csm(csc_2, dimension, nnz_2);
    convert_to_csc(dimension, nnz_2, matrix_2, csc_2);
    init_csm(csr_2, dimension, nnz_2);
    convert_to_csr(dimension, nnz_2, matrix_2, csr_2);

    init_dcsm(dcsc_1, dimension, nnz_1);
    convert_to_dcsc(dimension, nnz_1, csc_1, dcsc_1);
    init_dcsm(dcsr_1, dimension, nnz_1);
    convert_to_dcsr(dimension, nnz_1, csr_1, dcsr_1);
    init_dcsm(dcsc_2, dimension, nnz_2);
    convert_to_dcsc(dimension, nnz_2, csc_2, dcsc_2);
    init_dcsm(dcsr_2, dimension, nnz_2);
    convert_to_dcsc(dimension, nnz_2, csr_2, dcsr_2);


    struct Intersection *intersections = safe_calloc(dimension, sizeof(struct Intersection));
    int intersections_count = find_intersections(intersections, dcsc_1, dcsr_2, nnz_1, nnz_2);
    print_intersections(intersections, intersections_count); // TEST
    print_matrix(dimension, matrix_1); //TEST
    print_csm_matrix(dimension, nnz_1, csc_1->J, csc_1->I, csc_1->NUM); // TEST
    print_matrix(dimension, matrix_2); //TEST
    print_csm_matrix(dimension, nnz_2, csr_2->J, csr_2->I, csr_2->NUM); // TEST

    printf("INTERSECTIONS COUNT %d\n", intersections_count);
    struct HeapNode *temp_heap = (struct HeapNode*)safe_calloc(intersections_count, sizeof(struct HeapNode));
    struct HeapNode *result_heap = (struct HeapNode*)safe_calloc(min(nnz_1, dimension) * min(nnz_2, dimension), sizeof(struct HeapNode));
    int temp_heap_size = 0;
    int result_heap_size = 0;

    for(int i = 0; i < intersections_count; i++){
        mult_intersection_1(dcsc_1, dcsc_2, intersections, temp_heap, &temp_heap_size, &i);
        increment_list(intersections, i);
    }
    int emptyVar; // TEST
    print_heap(temp_heap, temp_heap_size);
    while(intersections_empty(intersections, intersections_count) == false){
        print_intersections(intersections, intersections_count);
        printf("In main: ");
        print_heap(temp_heap, temp_heap_size);
        struct HeapNode *current_node = pop_head(temp_heap, &temp_heap_size);
        printf("Extraxted node index = %d, value = %f\n", current_node->index, current_node->value);
        print_heap(temp_heap, temp_heap_size);
        if(result_heap_size == 0 || (result_heap[0].key[0] != current_node->key[0] || result_heap[0].key[1] != current_node->key[1])){
            insert(result_heap, &result_heap_size, current_node, &current_node->index);
        }else{
            result_heap[0].value += current_node->value;
        }
        if(intersections[current_node->index].empty == false){
            printf("TI TEST 1: ti = %d\n", current_node->index);
            mult_intersection_1(dcsc_1, dcsc_2, intersections, temp_heap, &temp_heap_size, &current_node->index);
            printf("TI TEST 4: ti = %d\n", current_node->index);
            increment_list(intersections, current_node->index);
        }
        print_heap(result_heap, result_heap_size);
        printf("TI TEST 5: ti = %d\n", current_node->index);
        scanf("%d", &emptyVar); // TEST
    }
    for(int i = 0; i < result_heap_size; i++){
        printf("row: %d; column: %d; value: %f\n", result_heap[i].key[0], result_heap[i].key[1], result_heap[i].value);
    }

    free(result_heap); result_heap = NULL;
    free(temp_heap); temp_heap = NULL;

    delete_dcsm(dcsr_2);
    delete_dcsm(dcsc_2);
    delete_dcsm(dcsr_1);
    delete_dcsm(dcsc_1);
    delete_csm(csr_2);
    delete_csm(csc_2);
    delete_csm(csr_1);
    delete_csm(csc_1);
}