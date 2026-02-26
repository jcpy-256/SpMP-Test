/**
 * 矩阵的转换操作
 * 包括矩阵的转置、基本的格式转换
*/

#ifndef TRANSPORSITION_H
#define TRANSPORSITION_H

#include<stdio.h>
#include<stdlib.h>

#include"matrix.h"
#include "dataStruc.h"
#include "utils.h"
#include "btools.h"


void pre_convert_to_lower_triangular(COO_MTX *coo)
{
    int lower_nnz = 0;
    int nodes = coo->nrow;
    for (int i = 0; i < coo->nnz; i++)
    {
        if (coo->col_idx[i] < coo->row_idx[i])
        {
            lower_nnz++;
        }
    }
    // printf("non diagonal lower triangular nnz: %d\n", lower_nnz);
    lower_nnz += coo->nrow;
    printf(" lower triangular nnz: %d\n", lower_nnz);

    int *col_idxices = (int *)g_malloc(sizeof(int) * lower_nnz);
    int *row_idxices = (int *)g_malloc(sizeof(int) * lower_nnz);
    double *data = (double *)g_malloc(sizeof(double) * lower_nnz);
    int index = 0;
    int col = 0;
    for (int i = 0; i < coo->nnz; i++)
    {
        if (col == coo->col_idx[i])
        {
            col_idxices[index] = col;
            row_idxices[index] = col;
            data[index] = 1.0;
            index++;
            col++;
        }
        if (coo->col_idx[i] < coo->row_idx[i])
        {
            col_idxices[index] = coo->col_idx[i];
            row_idxices[index] = coo->row_idx[i];
            data[index] = coo->val[i];
            // printf("(%d, %d) %lf\n",row_idxices[index],col_idxices[index],data[index]);
            index++;
        }
    }
    // printf("index:%d\n", index);
    // printf("col:%d, col->ncol:%d\n", col, coo->ncol);
    if (col < coo->ncol)
    {
        for (col; col < coo->ncol; col++)
        {
            col_idxices[index] = col;
            row_idxices[index] = col;
            data[index] = 1.0;
            index++;
            // printf("index:%d\n",index);
            // col++;
        }
    }
    // for (int i = 0; i < coo->nrow; i++)
    // {
    //     col_idxices[index] = i;
    //     row_idxices[index] = i;
    //     data[index] = 1.0;
    //     index++;
    // }
    // printf("index:%d\n", index);
    /* free old memory */
    free(coo->col_idx);
    free(coo->row_idx);
    free(coo->val);

    coo->col_idx = col_idxices;
    coo->row_idx = row_idxices;
    coo->val = data;
    coo->nnz = lower_nnz;
}

// void add_diag_element(CSR_MTX *csr)
// {
//     int nrow = csr->nrow;
//     int ncol = csr->ncol;
//     int nnz = csr->row_ptr[nrow];
//     int *row_ptr = csr->row_ptr;
//     int *col_idx = csr->col_idx;
//     int *new_row_ptr = (int *) g_calloc(nrow + 1, sizeof(int));
//     int *new_col_idx = (int *) g_calloc(nnz + nrow, sizeof(int));
//     double *new_val = (double*) g_malloc((nnz + nrow) * sizeof(double));
//     for (int ir = 0; ir < nrow; ir++)
//     {
//         bool flag = false;
//         for (int  beg = row_ptr[ir]; beg < row_ptr[ir]; ir++)
//         {
//             if (col_idx[beg] = ir) // 对角位置存在数据
//             {
//                 flag = true;
//                 continue;
//             }
            
//         }
//         if(!flag) // 对角位置不存在数据
//         {

//         }
//     }
// }



void transposition_COO_to_CSR(const COO_MTX *coo, CSR_MTX *csr);

void transposition_CSC_to_CSR(const int nnz, const CSC_MTX* csc, CSR_MTX* csr);


/**
 * @brief: matrix format transposition COO to CSR
 */
void transposition_COO_to_CSR(const COO_MTX *coo, CSR_MTX *csr)
{
    int nodes = coo->nrow;
    int nnz = coo->nnz;
    int *row_ptr = (int *)g_calloc(nodes + 1, sizeof(int));
    int *col_idx = (int *)g_calloc(nnz, sizeof(int));
    double *data = (double *)g_malloc(nnz * sizeof(double));

    // histogram in row pointer
    for (int i = 0; i < nnz; i++)
    {
        row_ptr[coo->row_idx[i] + 1]++;
    }

    prefix_sum(row_ptr, nodes + 1);
    int *csrRowIncr = (int *)g_calloc(nodes + 1, sizeof(int));
    memcpy(csrRowIncr, row_ptr, sizeof(int) * (nodes + 1));

    for (int i = 0; i < nnz; i++)
    {
        int row = coo->row_idx[i];
        col_idx[csrRowIncr[row]] = coo->col_idx[i];
        data[csrRowIncr[row]] = coo->val[i];
        csrRowIncr[row]++;
    }

    csr->col_idx = col_idx;
    csr->row_ptr = row_ptr;
    csr->val = data;
    csr->ncol = coo->ncol;
    csr->nrow = coo->nrow;

    #pragma omp parallel for
    for (int i=0; i < nnz; i++){
      quick_sort_key_val_pair(col_idx + row_ptr[i], data + row_ptr[i], row_ptr[i+1] - row_ptr[i] - 1);
    //   assert(is_sorted(col_idx + row_ptr[i], col_idx + row_ptr[i+1]));
    }
}

void transposition_CSC_to_CSR(const int nnz, const CSC_MTX* csc, CSR_MTX* csr)
{
    int nodes = csc->nrow;
    int* row_ptr = (int*) g_calloc(nodes+1, sizeof(int));
    int* col_idx = (int*) g_calloc(nnz, sizeof(int));
    double* data = (double*) g_malloc(nnz* sizeof(double));

    //histogram in row pointer
    for(int i=0; i< nnz; i++)
    {
        row_ptr[csc->row_idx[i]+1]++;
    }
    prefix_sum(row_ptr, nodes + 1);
    int* csrRowIncr = (int*) g_calloc(nodes+1, sizeof(int));
    memcpy(csrRowIncr, row_ptr, sizeof(int) * (nodes+1));

    /* insert nnz to CSR*/
    for(int col = 0; col < csc->ncol; col++)
    {
        for(int j = csc->col_ptr[col]; j< csc->col_ptr[col+1]; j++)
        {
            int row = csc->row_idx[j];
            col_idx[csrRowIncr[row]] = col;
            data[csrRowIncr[row]] = csc->val[j];
            csrRowIncr[row]++;
        }
    }
}



void matrix_transposition_CSR(const CSR_MTX *csr, CSR_MTX *new_csr,const int nnz)
{
    
}


#endif