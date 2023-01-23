#include <stdlib.h>
#include <stdbool.h>

struct CSMatrix{
    int *J;
    int *I;
    double *NUM;
};

struct DCSMatrix{
    int *AUX;
    int *J;
    int *P;
    int *I;
    double *NUM;
};

// Compressed Sparse Column representation of the matrix
void convert_to_csc(
    int dimension,
    int nnz,
    double **matrix,
    struct CSMatrix *csc)
{
    int column_ptr = 0;
    for(int column_idx = 0; column_idx < dimension; column_idx++){
        (*csc).J[column_idx] = column_ptr;
        for(int row_idx = 0; row_idx < dimension; row_idx++){
            if(compare_double(matrix[column_idx][row_idx], 0.0) == true) continue;
            (*csc).I[column_ptr] = row_idx;
            (*csc).NUM[column_ptr] = matrix[column_idx][row_idx];
            column_ptr++;
        }
    }
    (*csc).J[dimension] = column_ptr;
}

// Compressed Sparse Row representation of the matrix
void convert_to_csr(
    int dimension,
    int nnz,
    double **matrix,
    struct CSMatrix *csr)
{
    int row_ptr = 0;
    for(int row_idx = 0; row_idx < dimension; row_idx++){
        (*csr).J[row_idx] = row_ptr;
        for(int column_idx = 0; column_idx < dimension; column_idx++){
            if(compare_double(matrix[column_idx][row_idx], 0.0) == true) continue;
            (*csr).I[row_ptr] = column_idx;
            (*csr).NUM[row_ptr] = matrix[column_idx][row_idx];
            row_ptr++;
        }
    }
    (*csr).J[dimension] = row_ptr;
}

// Double Compressed Sparse Column representation of the matrix
void convert_to_dcsc(
    int dimension, int nnz,
    struct CSMatrix *csc,
    struct DCSMatrix *dcsc
    //int IC[dimension],
    //int CP[dimension + 1],
    //int JC[nnz + 1],
    //int AUX[nnz + 1]
){
    int nzc = (*csc).J[dimension];
    int cf = (dimension + 1) / nzc;

    int p_count = 0;
    int aux_count = 0;
    int column_idx = 0;
    bool aux_updated = false;
    while(column_idx < dimension)
    {
        aux_updated = false;
        for(int chunk_idx = 0; chunk_idx < cf; chunk_idx++)
        {
            column_idx++;
            int diff = (*csc).J[column_idx] - (*csc).J[column_idx - 1];
            if(diff == 0) continue;
            (*dcsc).J[p_count] = column_idx - 1;
            (*dcsc).P[p_count] = (*csc).J[column_idx - 1];
            if(aux_updated == false)
            {
                aux_updated = true;
                (*dcsc).AUX[aux_count] = p_count;
                aux_count++;
            }
            p_count++;
            if(column_idx == dimension) break;
        }
    }
    (*dcsc).P[p_count] = (*csc).J[dimension];
    if(aux_updated == false){
        (*dcsc).AUX[aux_count] = p_count;
        aux_count = aux_count + 1;
    }
    (*dcsc).AUX[aux_count] = (*dcsc).AUX[aux_count - 1];
    (*dcsc).I = (*csc).I;
    (*dcsc).NUM = (*csc).NUM;
}

// Double Compressed Sparse Row representation of the matrix
void convert_to_dcsr(
    int dimension, int nnz,
    struct CSMatrix *csr,
    struct DCSMatrix *dcsr
    //int IR[dimension],
    //int RP[dimension + 1],
    //int JR[nnz + 1],
    //int AUX[nnz + 1]
)
{
    int nzr = (*csr).J[dimension];
    int cf = (dimension + 1) / nzr + 1; // TODO: convert to int

    int p_count = 0;
    int aux_count = 0;
    int row_idx = 0;
    bool aux_updated = false;
    while(row_idx < dimension)
    {
        aux_updated = false;
        for(int chunk_idx = 0; chunk_idx < cf; chunk_idx++)
        {
            row_idx++;
            int diff = (*csr).J[row_idx] - (*csr).J[row_idx - 1];
            if(diff == 0) continue;
            (*dcsr).J[p_count] = row_idx;
            (*dcsr).P[p_count] = (*csr).J[row_idx];
            if(aux_updated == true)
            {
                aux_updated = true;
                (*dcsr).AUX[aux_count] = p_count;
                aux_count++;
            }
            p_count++;
            if(row_idx == dimension) break;
        }
    }
    (*dcsr).P[p_count] = (*csr).J[dimension];
    if(aux_updated == false){
        (*dcsr).AUX[aux_count] = p_count;
        p_count++;
        aux_count++;
    }
    (*dcsr).AUX[aux_count] = p_count - 1;
    copy_int_vector((*dcsr).I, (*csr).I, nnz);
    copy_double_vector((*dcsr).NUM, (*csr).NUM, nnz);
}


void init_csm(struct CSMatrix *csm, int dimension, int nnz){
    if (csm == NULL){
        printf("Fatal: failed to allocate %zu bytes.\n", sizeof(struct CSMatrix));
        abort();
    }
    (*csm).J = safe_calloc_int(dimension + 1);
    (*csm).I = safe_calloc_int(nnz);
    (*csm).NUM = safe_calloc_double(nnz);
}

void init_dcsm(struct DCSMatrix *dcsm, int dimension, int nnz){
    if (dcsm == NULL){
        printf("Fatal: failed to allocate %zu bytes.\n", sizeof(struct DCSMatrix));
        abort();
    }
    (*dcsm).AUX = safe_calloc_int(nnz + 1);
    (*dcsm).I = safe_calloc_int(dimension + 1);
    (*dcsm).J = safe_calloc_int(dimension);
    (*dcsm).NUM = safe_calloc_double(nnz);
    (*dcsm).P = safe_calloc_int(dimension + 1);
}

void delete_csm(struct CSMatrix *csm){
    free((*csm).I); (*csm).I = NULL;
    free((*csm).J); (*csm).J = NULL;
    free((*csm).NUM); (*csm).NUM = NULL;
    free(csm); csm = NULL;
}

void delete_dcsm(struct DCSMatrix *dcsm){
    free((*dcsm).AUX); (*dcsm).AUX = NULL;
    free((*dcsm).I); (*dcsm).I = NULL;
    free((*dcsm).J); (*dcsm).J = NULL;
    free((*dcsm).NUM); (*dcsm).NUM = NULL;
    free((*dcsm).P); (*dcsm).P = NULL;
    free(dcsm); dcsm = NULL;
}

struct HeapNode{
    int key[2];
    double value;
    int index;
};

void copy_heap_node(struct HeapNode *destination, struct HeapNode *source){
    destination->index = source->index;
    destination->key[0] = source->key[0];
    destination->key[1] = source->key[1];
    destination->value = source->value;
}

void swap(struct HeapNode *a, struct HeapNode *b, int *dti){
    printf("TI TEST swap.1: ti = %d\n", *dti);
    struct HeapNode temp = *b;
    printf("TI TEST swap.2: ti = %d\n", *dti);
    *b = *a;
    printf("TI TEST swap.3: ti = %d\n", *dti);
    *a = temp;
    printf("TI TEST swap.4: ti = %d\n", *dti);
}

void min_heapify(struct HeapNode *heap, int heap_size, int current, int *dti){
    //printf("TI TEST min_heapify.1: ti = %d\n", *dti);

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
        printf("TI TEST min_heapify.1.1: ti = %d, heap_size = %d\n", *dti, heap_size);
        printf("current = %d, smallest = %d\n", smallest, current);
        swap(&(heap[current]), &(heap[smallest]), dti);
        printf("TI TEST min_heapify.1.2: ti = %d\n", *dti);
        min_heapify(heap, heap_size, smallest, dti);
    }
    //printf("TI TEST min_heapify.2: ti = %d\n", *dti);
}

void insert(struct HeapNode *heap, int *heap_size, struct HeapNode *new_node, int *ti){
    printf("INSERTING key [%d, %d]\n", new_node->key[0], new_node->key[1]);
    printf("TI TEST 2.2.1: ti = %d\n", *ti);
    if(*heap_size == 0){
        copy_heap_node(&heap[0], new_node);
        (*heap_size)++;
        return;
    }
    copy_heap_node(&heap[*heap_size], new_node);
    (*heap_size)++;

    printf("TI TEST 2.2.2: ti = %d\n", *ti);
    for(int hi = (*heap_size) / 2 - 1; hi>= 0; hi--){
        min_heapify(heap, *heap_size, hi, ti);
        printf("TI TEST 2.2.2.%d: ti = %d\n", hi, *ti);
    }
    //printf("TI TEST 2.2.3: ti = %d\n", *ti);
}

void print_heap(struct HeapNode *heap, int heap_size){
    for(int i = 0; i < heap_size; i++){
        printf(
            "index = %d, key = [%d, %d], value = %f\n",
            heap[i].index, heap[i].key[0], heap[i].key[1], heap[i].value
        );
    }
    printf("\n");
}

struct HeapNode *pop_head(struct HeapNode *heap, int *heap_size){
    printf("In pop Head 1: \n");
    print_heap(heap, *heap_size);
    swap(&heap[0], &heap[(*heap_size) - 1], heap_size);
    printf("In pop Head 2: \n");
    print_heap(heap, *heap_size);
    (*heap_size)--;
    for(int i = (*heap_size) / 2 - 1; i >= 0; i--){
        min_heapify(heap, *heap_size, i, &i);
    }
    return &heap[(*heap_size)];
}