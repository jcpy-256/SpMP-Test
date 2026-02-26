/**
 * 2025-09-08 by yansen
 */
#include <iostream>
#include <cassert>
#include <cstring>
#include <climits>
#include <cfloat>

#include <omp.h>

#ifdef MKL
#include <mkl.h>
#endif

#include "../LevelSchedule.hpp"
#include "../synk/barrier.hpp"

#include "test.hpp"

#include "btools.h"
#include "sparseMtxOp.h"
#include "../COO.hpp"

using namespace SpMP;

#define CHECK_POINTER(ptr)                                                     \
    do                                                                         \
    {                                                                          \
        if ((ptr) == NULL)                                                     \
        {                                                                      \
            fprintf(stderr, "Error: Null pointer detected in %s at line %d\n", \
                    __FILE__, __LINE__);                                       \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0);

// parallel ilu0
void ilu0(double *lu, const CSR &A, const LevelSchedule &schedule)
{
    int base = A.getBase();

    const int *rowptr = A.rowptr - base;
    const int *colidx = A.colidx - base;
    const int *diagptr = A.diagptr - base;
    const double *values = A.values - base;

    lu -= base;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        // #pragma omp for
        //         for (int i = base; i < A.getNnz() + base; i++)
        //         {
        //             lu[i] = values[i];
        //         }

        const int ntasks = schedule.ntasks;
        const short *nparents = schedule.nparentsForward;
        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        const int *perm = schedule.threadContToOrigPerm;

        int nBegin, nEnd;
        getSimpleThreadPartition(&nBegin, &nEnd, ntasks);

        volatile int *taskFinished = schedule.taskFinished;
        int **parents = schedule.parentsForward;

        memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

        synk::Barrier::getInstance()->wait(tid);

        for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task)
        {
            SPMP_LEVEL_SCHEDULE_WAIT;

            for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i)
            {
                int row = perm[i] + base;

                for (int j = rowptr[row]; j < diagptr[row]; ++j)
                {
                    int c = colidx[j];
                    double tmp = lu[j] /= lu[diagptr[c]];

                    int k1 = j + 1, k2 = diagptr[c] + 1;

                    while (k1 < rowptr[row + 1] && k2 < rowptr[c + 1])
                    {
                        if (colidx[k1] < colidx[k2])
                            ++k1;
                        else if (colidx[k1] > colidx[k2])
                            ++k2;
                        else
                        {
                            lu[k1] -= tmp * lu[k2];
                            ++k1;
                            ++k2;
                        }
                    }
                }
            } // for each row

            SPMP_LEVEL_SCHEDULE_NOTIFY;
        } // for each level
    } // omp parallel
}

void ilu0csr_uplooking_ref(CSR *A, double *lu)
{

#define ILU0_INITIAL
    assert(A->diagptr);
    const int nnz = A->getNnz();
    const int m = A->m;
    const int n = A->n;
    const int *rowptr = A->rowptr;
    const int *colidx = A->colidx;
    const int *diagptr = A->diagptr;
    double *val = A->values;
    // if (m != n)
    // {
    //     throw std::runtime_error("this matrix is not a square matrix");
    // }
    for (int i = 0; i < m; i++)
    {
        for (int j = rowptr[i]; j < diagptr[i]; ++j) // processing only lower triangular part
        {
            int col = colidx[j];
            double pivot = lu[diagptr[col]];
            assert(fabs(pivot) > 1.0e-16);
            if (fabs(pivot) < 1.0e-16)
            {
                throw std::runtime_error("Zero pivot encountered at row " + std::to_string(col));
            }
            double tmp = lu[j] /= pivot;

            int k1 = j + 1, k2 = diagptr[col] + 1;
            while (k1 < rowptr[i + 1] && k2 < rowptr[col + 1])
            {
                if (colidx[k1] < colidx[k2])
                    k1++;
                else if (colidx[k1] > colidx[k2])
                    k2++;
                else
                {
                    lu[k1] -= tmp * lu[k2];
                    ++k1;
                    ++k2;
                }
            } // two pointer scan and update vector
        } // for: from rowptr[row] to diagptr[row], processing only lower triangular part
    } // for: traverse all row
}

CSR *forceLowerSymmmetric(const CSR *csrA)
{

    const int m = csrA->m;
    const int n = csrA->n;
    const int *csrColIdx = csrA->colidx;
    const int *csrRowPtr = csrA->rowptr;
    const double *csrVal = csrA->values;
    if (m != n)
    {
        printf("the matrix is no a square matrix, and the force symmetric couldn't proceed!\n");
        exit(1);
    }
    int nnz = csrA->getNnz();

    int *csrRowCounter = MALLOC(int, m + 1);
    // std::copy(csrA->rowptr, csrA->rowptr + m + 1, csrRowCounter);
    memcpy(csrRowCounter, csrA->rowptr, (m + 1) * sizeof(int));

    for (int i = 0; i < m; i++)
    {
        csrRowCounter[i] = csrRowCounter[i + 1] - csrRowCounter[i];
    }

    double maxAbsoluteValue = __DBL_MIN__;

    for (int i = 0; i < m; i++)
    {
        for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
        {
            if (i != csrColIdx[j])
            {
                csrRowCounter[csrColIdx[j]]++;
            }
        }
    }

    exclusive_scan(csrRowCounter, m + 1);
    // prefix_sum(csrRowCounter, m+1);
    // assert(csrRowCounter[m] == )

    int nnz_tmp = csrRowCounter[m]; // 0-based
    assert(nnz_tmp == 2 * nnz - m);
    int *csrRowPtr_alias = (int *)malloc((m + 1) * sizeof(int));
    CHECK_POINTER(csrRowPtr_alias);
    int *csrColIdx_alias = (int *)malloc(nnz_tmp * sizeof(int));
    CHECK_POINTER(csrColIdx_alias);
    double *csrVal_alias = (double *)malloc(nnz_tmp * sizeof(double));
    CHECK_POINTER(csrVal_alias);

    memcpy(csrRowPtr_alias, csrRowCounter, (m + 1) * sizeof(int));
    memset(csrRowCounter, 0, (m + 1) * sizeof(int));

    for (int i = 0; i < m; i++)
    {
        for (int j = csrRowPtr[i]; j < csrRowPtr[i + 1]; j++)
        {
            if (i != csrColIdx[j]) // not a diagonal element
            {
                // 原始值
                int row = i;
                int col = csrColIdx[j];
                assert(row > col);
                // val =
                int offset = csrRowPtr_alias[row] + csrRowCounter[row];
                csrColIdx_alias[offset] = col;
                csrVal_alias[offset] = csrVal[j];
                csrRowCounter[row]++;
                // 对称值

                row = csrColIdx[j];
                col = i;
                assert(col > row);
                offset = csrRowPtr_alias[row] + csrRowCounter[row];
                csrColIdx_alias[offset] = col;
                csrVal_alias[offset] = csrVal[j];
                csrRowCounter[row]++;
            }
            else // diagonal element
            {
                // assert(row == col);
                int offset = csrRowPtr_alias[i] + csrRowCounter[i];
                csrColIdx_alias[offset] = csrColIdx[j];
                csrVal_alias[offset] = csrVal[j];
                // csrVal_alias[offset] = maxAbsoluteValue;
                csrRowCounter[i]++;
            }
        }
    }

#pragma omp parallel for
    for (int i = 0; i < m; i++)
    {
        // qsort(csrColIdx_alias, csrVal_alias, csrRowPtr_alias[i], csrRowPtr_alias[i + 1] - 1);
        quicksort(csrColIdx_alias, csrVal_alias, csrRowPtr_alias[i], csrRowPtr_alias[i + 1] - 1);
    }

    CSR *csrTmp = new CSR();
    csrTmp->m = m;
    csrTmp->n = n;
    csrTmp->colidx = csrColIdx_alias;
    csrTmp->rowptr = csrRowPtr_alias;
    csrTmp->values = csrVal_alias;

    csrTmp->diag = MALLOC(double, m);
    csrTmp->idiag = MALLOC(double, m);
    printf("orig m:%d, n:%d, nnz:%d\n", m, n, nnz);

    fflush(stdout);
    csrTmp->constructDiagPtr();
    // csrTmp->computeInverseDiag()

    CSR *csrSym = new CSR(*csrTmp);
    delete csrTmp;

    return csrSym;
}

// serial ilu0
void ilu0_ref(double *lu, const CSR &A)
{
    int base = A.getBase();

    const int *rowptr = A.rowptr - base;
    const int *colidx = A.colidx - base;
    const int *diagptr = A.diagptr - base;
    const double *values = A.values - base;

    lu -= base;

    // #pragma omp for
    //     for (int i = base; i < A.getNnz() + base; i++)
    //     {
    //         lu[i] = values[i];
    //     }

    for (int i = 0; i < A.m; ++i)
    {
        for (int j = rowptr[i]; j < diagptr[i]; ++j)
        {
            int c = colidx[j];
            double tmp = lu[j] /= lu[diagptr[c]];

            int k1 = j + 1, k2 = diagptr[c] + 1;

            while (k1 < rowptr[i + 1] && k2 < rowptr[c + 1])
            {
                if (colidx[k1] < colidx[k2])
                    ++k1;
                else if (colidx[k1] > colidx[k2])
                    ++k2;
                else
                {
                    lu[k1] -= tmp * lu[k2];
                    ++k1;
                    ++k2;
                }
            }
        }
    } // for each row
}

int read_matrix_list(char ***matrix_list, const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        printf("file %s is not exist!!\n", filename);
        exit(2);
    }
    char matrix_name[256];
    int total = 0;
    int max_len = 256;

    while (!feof(fp))
    {
        fscanf(fp, "%s\n", matrix_name);
        memset(matrix_name, 0, 256 * sizeof(char));
        // max_len = max_len > strlen(matrix_name) ? max_len : strlen(matrix_name);
        total++;
    }

    char **matrix_temp = (char **)g_malloc(total * sizeof(char *));

    for (int i = 0; i < total; i++)
    {
        matrix_temp[i] = (char *)g_malloc((max_len + 1) * sizeof(char));
    }
    // printf("read end\n");
    rewind(fp); // 重新将文件指针回到文件的开始
    int idx = 0;
    while (!feof(fp) && idx < total)
    {
        fscanf(fp, "%s\n", matrix_temp[idx++]);
        // printf("matrix: %s\n", matrix_temp[idx - 1]);
    }

    fclose(fp);

    *matrix_list = matrix_temp;
    return idx;
}

void write_res_to_csv(const char *filename,
                      const char *matrix_name,
                      const int test_turns,
                      const int nrow,
                      const int nnzA,
                      const int levels,
                      const double serial_time,
                      const double p2p_time)

{
    char header[4096] = {'\0'};
    strcpy(header, "matrix");
    strcat(header, ",test_timestamp");
    strcat(header, ",test_turns");

    strcat(header, ",nrow");
    strcat(header, ",nnzA");
    strcat(header, ",levels");
    strcat(header, ",serial_time");
    strcat(header, ",p2p_time");
    char *content = (char *)g_malloc(4048 * sizeof(char));
    memset(content, 0, 4048 * sizeof(char));
    char cache[256] = {'\0'};
    // strcat(header, ",test_timestamp");
    char time_str[64] = {'\0'};
    time_t now_time;
    struct tm *p;
    time(&now_time);
    p = gmtime(&now_time);
    sprintf(time_str, "%d-%02d-%02d %02d:%02d:%02d", 1900 + p->tm_year, 1 + p->tm_mon, p->tm_mday, 8 + p->tm_hour, p->tm_min, p->tm_sec);
    sprintf(content, "%s,%s,%ld,%ld,%ld,%ld,%lf,%lf",
            matrix_name, time_str, test_turns, nrow, nnzA, levels, serial_time, p2p_time);

    if (fileIsExist(filename) == 0)
    {
        wirte_item_to_csv("\0", content, filename);
    }
    else
    {
        wirte_item_to_csv(header, content, filename);
    }
    free(content);
}

/** 将下三角矩阵进行对称转为整体的矩阵，对角位置元素不变 one-based */
void make_full(COO_MTX *coo)
{
    int *rowidx = coo->row_idx;
    int *colidx = coo->col_idx;
    double *values = coo->val;
    int m = coo->nrow;
    int org_nnz = coo->nnz;
    int count = (coo->nnz - m) * 2 + m;

    int *rowidx_new = (int *)malloc(count * sizeof(int));
    int *colidx_new = (int *)malloc(count * sizeof(int));
    double *values_new = (double *)malloc(count * sizeof(double));
    CHECK_POINTER(rowidx_new)
    CHECK_POINTER(colidx_new)
    CHECK_POINTER(values_new)
    //    FREE
    memcpy(rowidx_new, rowidx, org_nnz * sizeof(int));
    memcpy(colidx_new, colidx, org_nnz * sizeof(int));
    memcpy(values_new, values, org_nnz * sizeof(double));

    for (int i = 0; i < coo->nnz; i++)
    {
        int row = rowidx[i];
        int col = colidx[i];
        if (row > col)
        {
            rowidx_new[org_nnz] = col;
            colidx_new[org_nnz] = row;
            values_new[org_nnz] = values[i];
            org_nnz++;
        }
    }
    if (org_nnz != count)
    {
        printf("symmetric opertaion failure in %s of file %s\n", __LINE__, __FILE__);
        exit(2);
    }
    coo->col_idx = colidx_new;
    coo->row_idx = rowidx_new;
    coo->val = values_new;
    coo->nnz = count;

/* make one-based*/
#pragma omp parallel for
    for (int i = 0; i < count; i++)
    {
        coo->col_idx[i]++;
        coo->row_idx[i]++;
    }

    free(rowidx);
    free(colidx);
    free(values);
}

void ILU0_kernel(const char *matrix_name, const char *csv_filename, int num_test)
{
    printf("file path:%s\n", matrix_name);
    int m, n;
    int nnz;
    int isSymmetric;

    CSR *A = new CSR(matrix_name);
    m = A->m, n = A->n;
    nnz = A->getNnz();
    // printf("CSR m: %d, n: %d, nnz: %d\n", m, n, nnz);

    // constrcut lower matrix to get schedule
    int lower_nnz = 0;
    for (int i = 0; i < m; i++)
    {
        lower_nnz += A->diagptr[i] - A->rowptr[i] + 1;
    }
    printf("lower triangular part has %d nnz\n", lower_nnz);
    fflush(stdout);

    CSR * csrALower = new CSR(m,n,lower_nnz);
    csrALower->rowptr[0] = 0;
    for (int i = 0; i < m; i++)
    {
        csrALower->rowptr[i + 1] = csrALower->rowptr[i] + A->diagptr[i] - A->rowptr[i] + 1;
        std::copy(A->colidx + A->rowptr[i], A->colidx + A->diagptr[i] + 1, csrALower->colidx + csrALower->rowptr[i]);
    }

    assert(csrALower->rowptr[m] == lower_nnz);
    assert(csrALower->colidx[lower_nnz-1] == m - 1);

    // // make_full()
    CSR* csrALowerSym = forceLowerSymmmetric(csrALower);

    // ************************** Serial ILU0 Test ******************************
    std::vector<double> timeVec(num_test, 0.0);
    double *lu_serial = MALLOC(double, nnz);

    double a = getCurrentTimeMilli();
    for (int i = 0; i < num_test; i++)
    {
        // copy value to lu_serial
        copyVector(lu_serial, A->values, nnz);
        double t1 = getCurrentTimeMilli();
        // ilu0_ref(lu_serial, *A);
        ilu0csr_uplooking_ref(A,lu_serial);
        timeVec[i] = getCurrentTimeMilli() - t1;
    }

    printf("serial time using is %lf\n", getCurrentTimeMilli() - a);
    fflush(stdout);
    // get median
    std::sort(timeVec.begin(), timeVec.end());
    double serial_time = num_test % 2 == 1 ? timeVec[(num_test) / 2] : (timeVec[(num_test + 1) / 2] + timeVec[(num_test) / 2]) / 2;

    // ************************** p2p ILU0 Test ******************************
    LevelSchedule p2pScheduleWithTransitiveReduction;
    int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
    bool wasSymmetric = getSymmetricNnzPattern(csrALowerSym, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);
    // printf("issymmetric: %s\n", isSymmetric ? "true" : "false");
    // p2pScheduleWithTransitiveReduction.constructTaskGraph(*A);
    if (wasSymmetric)
    {
        p2pScheduleWithTransitiveReduction.constructTaskGraph(*csrALowerSym);
    }
    else
    {
        // p2pScheduleWithTransitiveReduction.constructTaskGraph(
        //     ASym->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
        //     PrefixSumCostFunction(symRowPtr));
    }

    FREE(symRowPtr);
    FREE(symColIdx);
    FREE(symDiagPtr);
    FREE(symExtPtr);
    int nlevels = 0;
    printf("parallelism = %g\n", (double)A->m / (p2pScheduleWithTransitiveReduction.levIndices.size() - 1));
    nlevels = p2pScheduleWithTransitiveReduction.levIndices.size() - 1;
    double *lu_p2p = MALLOC(double, A->getNnz());
    for (int i = 0; i < num_test; i++)
    {
        // copy value to lu_serial
        copyVector(lu_p2p, A->values, nnz);
        double t1 = getCurrentTimeMilli();
        ilu0(lu_p2p, *A, p2pScheduleWithTransitiveReduction);
        timeVec[i] = getCurrentTimeMilli() - t1;
        // printf("turn: %d\n", i);
    }
    // get median
    std::sort(timeVec.begin(), timeVec.end());
    double p2p_time = num_test % 2 == 1 ? timeVec[(num_test) / 2] : (timeVec[(num_test + 1) / 2] + timeVec[(num_test) / 2]) / 2;

    write_res_to_csv(csv_filename, matrix_name,num_test, A->m, A->getNnz(), nlevels, serial_time, p2p_time);
    // write_res_to_csv(csv_filename, matrix_name, A->m, A->getNnz(), 0, serial_time, p2p_time);
    printf("serial time: %lf, p2p time: %lf\n", serial_time, p2p_time);
    bool check = correctnessCheck(lu_p2p, lu_serial, A->getNnz(), 1.0e-10);
    if (check)
    {
        printf("ILU check pass!\n");
    }
    else
    {
        printf("ILU check not pass!\n");
    }
    FREE(lu_serial);
    FREE(lu_p2p);
    delete A;
}

void run( const char *matrix_name,
         const char *csv_res_file, int num_test)
{
    printf("###################### %s #######################\n", matrix_name);
    ILU0_kernel(matrix_name, csv_res_file, num_test);
}

int main(int argc, char *argv[])
{
    int nthreads = atoi(argv[3]);
    int num_test = atoi(argv[4]);
    // printf("nthread: %d\n", nthreads);
    omp_set_num_threads(nthreads);

    char matrix_name[1024] = {'\0'};
    strcpy(matrix_name, argv[1]);
    // printf("matrix file: %s\n", matrix_name);
    char csv_res_file[1024] = "\0";
    strcpy(csv_res_file, argv[2]);
    run(matrix_name, csv_res_file, num_test);
}