/**
 * @file    : sparse_matrix.h
 * @author  : theSparky Team
 * @version :
 *
 * Functions for Sparse Matrix Init.
 */
#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>

struct COO_Matrix
{
	// public:
	int nrow;
	int ncol;
	int nnz;

	int *row_idx;
	int *col_idx;
	double *val;

	// COO_Matrix();
	// COO_Matrix(int n, int m, int nnz, int* row_ind, int* col_ind, double* values);
	// COO_Matrix(const COO_Matrix& A);
	// ~COO_Matrix();
	// COO_Matrix& operator=(const COO_Matrix& A);

	// void Free();
};

struct CSR_Matrix
{
	// public:
	int nrow;
	int ncol;

	int *row_ptr;
	int *col_idx;
	double *val;
	double *diagonal; // for SymGS
	int nnz;

	// CSR_Matrix();
	// CSR_Matrix(int n, int m, int* row_ptr, int* col_ind, double* values, double* diagonal);
	// CSR_Matrix(const CSR_Matrix& A);
	// CSR_Matrix(const COO_Matrix& A);
	// ~CSR_Matrix();
	// CSR_Matrix& operator=(const CSR_Matrix& A);
	// CSR_Matrix& operator=(const COO_Matrix& A);

	// void Free();
};

struct CSC_Matrix
{
	// public:
	int nrow;
	int ncol;

	int *row_idx;
	int *col_ptr;
	double *val;
	int nnz;

	// CSC_Matrix();
	// CSC_Matrix(int n, int m, int* row_ind, int* col_ptr, double* values);
	// CSC_Matrix(const CSC_Matrix& A);
	// CSC_Matrix(const COO_Matrix& A);
	// ~CSC_Matrix();
	// CSC_Matrix& operator=(const CSC_Matrix& A);
	// CSC_Matrix& operator=(const COO_Matrix& A);

	// void Free();
};

struct ELL_Matrix
{
	// public:
	int nrow;
	int ncol;
	int nnz;
	int nonzeros_in_row;

	int *col_idx;
	double *vals;
	double *diagonal; // for SymGS

	// ELL_Matrix();
	// ELL_Matrix(int n, int m, int nnz, int nonzeros_in_row, int* col_ind, double* values, double* diagonal);
	// ELL_Matrix(const ELL_Matrix& A);
	// ELL_Matrix(const COO_Matrix& A);
	// ~ELL_Matrix();
	// ELL_Matrix& operator=(const ELL_Matrix& A);
	// ELL_Matrix& operator=(const COO_Matrix& A);

	// void Free();
};

struct BSR_Matrix
{
	int blockDim;		// block dimsension of Martix
	int nrowb;
	int ncolb;
	int nnzb;
	double *bsrVal;	// length == nnzb * blockDim * blockDim
	int *bsrRowPtr;	// length == nrowb + 1
	int* bsrColIdx;	// length == nnzb

};




#endif