#include <stdbool.h>
#include <string.h>
#include "application/common/common.h"

typedef struct{
    int *J;
    int *I;
    double *NUM;
} CSMatrix;

CSMatrix *init_csm(int dimension, int nnz, double **matrix, bool column_major){
    CSMatrix *csm = malloc(sizeof(CSMatrix));
    if (csm == NULL){
        printf("Fatal: failed to allocate %zu bytes for csm.\n", sizeof(CSMatrix));
        abort();
    }
    (*csm).J = safe_calloc_int(dimension + 1);
    (*csm).I = safe_calloc_int(nnz);
    (*csm).NUM = safe_calloc_double(nnz);

    int major_ptr = 0;
    for(int major_idx = 0; major_idx < dimension; major_idx++){
        (*csm).J[major_idx] = major_ptr;
        for(int minor_idx = 0; minor_idx < dimension; minor_idx++){
            double test_value = column_major==true ? matrix[major_idx][minor_idx]: matrix[minor_idx][major_idx];
            if(compare_double(test_value, 0.0) == true) continue;
            (*csm).I[major_ptr] = minor_idx;
            (*csm).NUM[major_ptr] = test_value;
            major_ptr++;
        }
    }
    (*csm).J[dimension] = major_ptr;

    return csm;
}

void print_csm_matrix(int dimension, int nnz, CSMatrix *csm){
    for(int i = 0; i < dimension; i ++){
        printf("%.5i ", csm->J[i]);
        for (int j = 0; j < csm->J[i + 1] - csm->J[i] - 1; j++){
            printf("      ");
        }
    }
    printf("%.5i\n", csm->J[dimension]);
    for(int i = 0; i < nnz; i ++){
        printf("%.5i ", csm->I[i]);
    }
    printf("\n");
    for(int i = 0; i < nnz; i ++){
        printf("%5.2f ", csm->NUM[i]);
    }
    printf("\n");
}

void delete_csm(CSMatrix *csm){
    free((*csm).J);
    free((*csm).I);
    free((*csm).NUM);
    free(csm); csm = NULL;
}

typedef struct{
    int *AUX;
    int *J;
    int *P;
    int *I;
    double *NUM;
} DCSMatrix;

/*DCSMatrix *init_dcsm(int dimension, int nnz, CSMatrix *csm){
    DCSMatrix *dcsm = malloc(sizeof(DCSMatrix));
    if (dcsm == NULL){
        printf("Fatal: failed to allocate %zu bytes.\n", sizeof(DCSMatrix));
        abort();
    }
    int nz = csm->J[dimension];

    dcsm->AUX = safe_calloc_int(nz + 2);
    dcsm->J = safe_calloc_int(nz);
    dcsm->P = safe_calloc_int(nz + 1);
    dcsm->I = safe_calloc_int(nnz);
    dcsm->NUM = safe_calloc_double(nnz);
    int cf = (dimension + 1) / nz;
    cf = cf == 0? 1: cf;
    printf("cf = %d\n", cf);
    int p_count = 0;
    int aux_count = 0;
    int major_idx = 0;
    bool aux_updated = false;
    while(major_idx < dimension)
    {
        aux_updated = false;
        for(int chunk_idx = 0; chunk_idx < cf; chunk_idx++)
        {
            major_idx++;
            int diff = csm->J[major_idx] - csm->J[major_idx - 1];
            if(diff == 0) continue;
            dcsm->J[p_count] = major_idx - 1;
            dcsm->P[p_count] = csm->J[major_idx - 1];
            if(aux_updated == false)
            {
                aux_updated = true;
                (*dcsm).AUX[aux_count] = p_count;
                aux_count++;
            }
            p_count++;
            if(chunk_idx - 1 == cf) aux_updated = false;
            if(major_idx == dimension) break;
        }
    }
    dcsm->P[p_count] = csm->J[dimension];
    if(aux_updated == false){
        printf("Aux updated = false\n");
        dcsm->AUX[aux_count] = p_count;
        aux_count++;
    }
    dcsm->AUX[aux_count] = dcsm->AUX[aux_count - 1];
    copy_int_vector(dcsm->I, csm->I, nnz);
    copy_double_vector(dcsm->NUM, csm->NUM, nnz);

    return dcsm;
}*/

DCSMatrix *init_dcsm(int dimension, int nnz, CSMatrix *csm){
    DCSMatrix *dcsm = malloc(sizeof(DCSMatrix));
    if (dcsm == NULL){
        printf("Fatal: failed to allocate %zu bytes for dcsm.\n", sizeof(DCSMatrix));
        abort();
    }
    
    int nz = csm->J[dimension];

    dcsm->AUX = safe_calloc_int(nz + 2);
    dcsm->J = safe_calloc_int(nz);
    dcsm->P = safe_calloc_int(nz + 1);
    dcsm->I = safe_calloc_int(nnz);
    dcsm->NUM = safe_calloc_double(nnz);

    if(nz == 0) return dcsm;
    int cf = (dimension + 1) / nz;
    cf = cf == 0? 1: cf;
    int p_counter = 0;
    int aux_counter = 0;
    bool aux_updated = false;
    int i;
    for(i = 0; i < dimension; i++){
        int diff = csm->J[i + 1] - csm->J[i];
        if(diff != 0){
            dcsm->P[p_counter] = csm->J[i];
            dcsm->J[p_counter] = i;
            if(aux_updated == false){
                aux_updated = true;
                dcsm->AUX[aux_counter] = p_counter;
                aux_counter++;
            }
            p_counter++;
        }
        if((i + 1) % cf == 0) aux_updated = false;
    }
    dcsm->P[p_counter] = csm->J[dimension];
    if((i + 1) % cf == 0) aux_updated = false;
    if(aux_updated == false){
        aux_updated = true;
        dcsm->AUX[aux_counter] = p_counter;
        aux_counter++;
    }
    dcsm->AUX[aux_counter] = dcsm->AUX[aux_counter - 1];

    copy_int_vector(dcsm->I, csm->I, nnz);
    copy_double_vector(dcsm->NUM, csm->NUM, nnz);

    return dcsm;
}

void print_dcsm_matrix(DCSMatrix *dcsm, int nnz){
    int aux_i = 0;
    printf("AUX ");
    do{
        printf("%d, ", dcsm->AUX[aux_i]);
        aux_i++;
    } while(dcsm->AUX[aux_i] != dcsm->AUX[aux_i - 1]);
    printf("%d\n", dcsm->AUX[aux_i]);
    printf("  J ");
    for(int i = 0; i < aux_i - 2; i++){
        printf("%d, ", dcsm->J[i]);
    }
    printf("%d\n", dcsm->J[aux_i - 2]);
    printf("  P ");
    for(int i = 0; i < aux_i - 1; i++){
        printf("%d, ", dcsm->P[i]);
    }
    printf("%d\n", dcsm->P[aux_i - 1]);
    printf("  I ");
    for(int i = 0; i < nnz - 1; i++){
        printf("%d, ", dcsm->I[i]);
    }
    printf("%d\n", dcsm->I[nnz - 1]);
    printf("NUM ");
    for(int i = 0; i < nnz - 1; i++){
        printf("%.3f, ", dcsm->NUM[i]);
    }
    printf("%.3f\n\n", dcsm->NUM[nnz - 1]);
}

void delete_dcsm(DCSMatrix *dcsm){
    free((*dcsm).AUX); (*dcsm).AUX = NULL;
    free((*dcsm).I); (*dcsm).I = NULL;
    free((*dcsm).J); (*dcsm).J = NULL;
    free((*dcsm).NUM); (*dcsm).NUM = NULL;
    free((*dcsm).P); (*dcsm).P = NULL;
    free(dcsm); dcsm = NULL;
}

typedef struct{
    int row;
    int column;
    double value;
    int index;
} HeapNode;

HeapNode *init_heap(int count){
    HeapNode *heap = malloc(count * sizeof(HeapNode));
    if (heap == NULL){
        printf("Fatal: failed to allocate %zu bytes for heap.\n", count * sizeof(DCSMatrix));
        abort();
    }
    return heap;
}

void extend_heap(HeapNode *heap, int *allocations){
    printf("Reallocating %d -> %d\n", *allocations, (int)((*allocations + 10) * 1.1));
    *allocations = (int)((*allocations + 10) * 1.1);
    printf("Reallocating\n");
    heap = realloc(heap, *allocations * sizeof(HeapNode));
    if (heap == NULL){
        printf("Fatal: failed to REallocate %zu bytes for heap.\n", *allocations * sizeof(HeapNode));
        abort();
    }
}

void swap_heap_nodes(HeapNode *a, HeapNode *b){
    HeapNode temp = *b;
    *b = *a;
    *a = temp;
}

void min_heapify(HeapNode *heap, int heap_size, int current){
    if(heap_size <= 1) return;
    int lc = 2 * current + 1;
    int rc = 2 * current + 2;
    int smallest = current;

    if(lc < heap_size && heap[lc].value < heap[smallest].value){
        smallest = lc;
    }
    if(rc < heap_size && heap[rc].value < heap[smallest].value){
        smallest = rc;
    }
    if(current != smallest){
        swap_heap_nodes(&heap[current], &heap[smallest]);
        min_heapify(heap, heap_size, smallest);
    }
}

void insert_to_heap(HeapNode *heap, HeapNode *new_node, int *heap_size, int *allocations){
    if(*heap_size == *allocations){
        printf("heap size = %d\n", *heap_size);
        extend_heap(heap, allocations);
    }
    heap[*heap_size].column = new_node->column;
    heap[*heap_size].row = new_node->row;
    heap[*heap_size].value = new_node->value;
    heap[*heap_size].index = new_node->index;
    (*heap_size)++;

    if(*heap_size == 1) return;
    for(int hi = (*heap_size) / 2 - 1; hi>= 0; hi--){
        min_heapify(heap, *heap_size, hi);
    }
}

HeapNode pop_heap_head(HeapNode *heap, int *heap_size){
    HeapNode head;
    head.column = heap[0].column;
    head.row = heap[0].row;
    head.index = heap[0].index;
    head.value = heap[0].value;

    swap_heap_nodes(&heap[0], &heap[*heap_size - 1]);
    (*heap_size)--;
    for(int i = (*heap_size) / 2 - 1; i >= 0; i--){
        min_heapify(heap, *heap_size, i);
    }
    return head;
}

void print_heap(HeapNode *heap, int heap_size){
    for(int i = 0; i < heap_size; i++){
        printf(
            "index = %d, key = [%d, %d], value = %f\n",
            heap[i].index, heap[i].row, heap[i].column, heap[i].value
        );
    }
    printf("\n");
}

void delete_heap(HeapNode *heap){
    free(heap); heap = NULL;
}

typedef struct{
    int index;
    int sizeA;
    int sizeB;
    int startA;
    int startB;
    int curA;
    int curB;
    bool empty;
} Intersection;

Intersection *init_intersections(int count, int *intersections_count, DCSMatrix *dcsc, DCSMatrix *dcsr){
    Intersection *intersections = malloc(count * sizeof(Intersection));
    if (intersections == NULL){
        printf("Fatal: failed to allocate %zu bytes for intersections.\n", count * sizeof(Intersection));
        abort();
    }

    int nzc, nzr;
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
            intersections[*intersections_count].index = dcsc->J[idx_1];
            intersections[*intersections_count].sizeA = dcsc->P[idx_1 + 1] - dcsc->P[idx_1];
            intersections[*intersections_count].sizeB = dcsr->P[idx_2 + 1] - dcsr->P[idx_2];
            intersections[*intersections_count].startA = dcsc->P[idx_1];
            intersections[*intersections_count].startB = dcsr->P[idx_2];
            intersections[*intersections_count].curA = dcsc->P[idx_1];
            intersections[*intersections_count].curB = dcsr->P[idx_2];
            intersections[*intersections_count].empty = false;
            idx_1++;
            idx_2++;
            (*intersections_count)++;
        }else if(dcsc->J[idx_1] < dcsr->J[idx_2]){
            idx_1++;
        }else{
            idx_2++;
        }
    }
    return intersections;
}

void print_intersections(Intersection *intersections, int intersections_count){
    int printed_intersections = 0;
    int i = -1;
    printf("i, szA, szB, stA, stB, crA, crB, empty\n");
    while(printed_intersections != intersections_count){
        i++;
        if(intersections[i].empty == true) continue;
        printed_intersections++;
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

void increment_intersection(Intersection *intersections, int ti, int *intersections_count){
    if(intersections[ti].curB - intersections[ti].startB < intersections[ti].sizeB - 1){
        intersections[ti].curB = intersections[ti].curB + 1;
    }else{
        intersections[ti].curB = intersections[ti].startB;
        if(intersections[ti].curA - intersections[ti].startA < intersections[ti].sizeA - 1){
            intersections[ti].curA = intersections[ti].curA + 1;
        }else{
            (*intersections_count)--;
            intersections[ti].empty = true;
        }
    }
}

void delete_intersections(Intersection *intersections){
    free(intersections); intersections = NULL;
}

void get_next_intersection_result(DCSMatrix *dcsc, DCSMatrix *dcsr, HeapNode *heap, int *heap_size, int *alloc_heap, Intersection *intersections, int ti){
    int ptr_a = intersections[ti].curA;
    int ptr_b = intersections[ti].curB;
    double product = dcsc->NUM[ptr_a] * dcsr->NUM[ptr_b];
    HeapNode new_node;
    new_node.row = dcsc->I[ptr_a]; // TODO
    new_node.column = dcsr->I[ptr_b]; // TODO is it really like that?
    new_node.index = ti;
    new_node.value = product;
    insert_to_heap(heap, &new_node, heap_size, alloc_heap);
}

void hypersparse_gemm_1(
    int dimension, int nnz_1, int nnz_2,
    double **matrix_1,
    double **matrix_2,
    struct algorithm_result *result
){
    strcpy((*result).name, __FUNCTION__);

    CSMatrix *csc_1 = init_csm(dimension, nnz_1, matrix_1, true);
    CSMatrix *csr_2 = init_csm(dimension, nnz_2, matrix_2, false);
    DCSMatrix *dcsc_1 = init_dcsm(dimension, nnz_1, csc_1);
    DCSMatrix *dcsr_2 = init_dcsm(dimension, nnz_2, csr_2);
    
    // Starting the algorithm
    clock_t t = clock();

    int intersections_count = 0;
    Intersection *intersections = init_intersections(dimension, &intersections_count, dcsc_1, dcsr_2);
    
    HeapNode *temp_heap = init_heap(intersections_count);
    int temp_heap_size = 0;
    int alloc_temp_heap = intersections_count;

    HeapNode *result_heap = init_heap(dimension * dimension); // TODO optimize
    int result_heap_size = 0;
    int alloc_result_heap = dimension * dimension;

    int temp = intersections_count;
    for(int ti = 0; ti < temp; ti++)
    {
        get_next_intersection_result(dcsc_1, dcsr_2, temp_heap, &temp_heap_size, &alloc_temp_heap, intersections, ti);
        increment_intersection(intersections, ti, &intersections_count);
    }
    while(temp_heap_size != 0){
        HeapNode next_node = pop_heap_head(temp_heap, &temp_heap_size);

        if(result_heap_size == 0) insert_to_heap(result_heap, &next_node, &result_heap_size, &alloc_result_heap);
        else if(result_heap[0].column != next_node.column || result_heap[0].row != next_node.row) insert_to_heap(result_heap, &next_node, &result_heap_size, &alloc_result_heap);
        else result_heap[0].value += next_node.value;

        int ti = next_node.index;
        if(intersections[ti].empty == false){
            get_next_intersection_result(dcsc_1, dcsr_2, temp_heap, &temp_heap_size, &alloc_temp_heap, intersections, ti);
            increment_intersection(intersections, ti, &intersections_count);
        }
    }

    (*result).clock_ticks = (int)(clock() - t);
    for(int i = 0; i < result_heap_size; i++){
        result->matrix[result_heap[i].column][result_heap[i].row] += result_heap[i].value;
    }

    delete_heap(result_heap);
    delete_heap(temp_heap);
    delete_intersections(intersections);
    delete_dcsm(dcsr_2);
    delete_dcsm(dcsc_1);
    delete_csm(csr_2);
    delete_csm(csc_1);
}
