/**
 * 2025年4月14日
 */

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

#include "sparseMtxOp.h"
#include "../COO.hpp"

using namespace SpMP;

/**
 * Reference sequential sparse triangular solver
 */
void forwardSolveRef(const CSR &A, double y[], const double b[])
{
    ADJUST_FOR_BASE;
    double *diag = A.diag;
    for (int i = base; i < A.m + base; ++i)
    {
        double sum = b[i];
        for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
        {
            sum -= values[j] * y[colidx[j]];
        }
        y[i] = sum / diag[i];
    } // for each row
}

void backwardSolveRef(const CSR &A, double y[], const double b[])
{
    ADJUST_FOR_BASE;

    for (int i = A.m - 1 + base; i >= base; --i)
    {
        double sum = b[i];
        for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
        {
            sum -= values[j] * y[colidx[j]];
        }
        y[i] = sum;
    } // for each row
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization
 */
void forwardSolveWithBarrier(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule,
    const int *perm)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task)
        {
            for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i)
            {
                int row = perm[i] + base;
                double sum = b[row];
                for (int j = rowptr[row]; j < rowptr[row + 1]; ++j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[row] = sum * idiag[row];
            } // for each row

            synk::Barrier::getInstance()->wait(tid);
        } // for each level
    } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization
 */
void backwardSolveWithBarrier(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule,
    const int *perm)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;
        // 这里的task是按照level组织的，task的数量大致为nlevel * nthread 每个task结束时进行一次全局的同步
        for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task)
        {
            for (int i = taskBoundaries[task + 1] - 1; i >= taskBoundaries[task]; --i)
            {
                int row = perm[i] + base;
                double sum = b[row];
                for (int j = rowptr[row + 1] - 1; j >= rowptr[row]; --j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[row] = sum;
            } // for each row
            synk::Barrier::getInstance()->wait(tid);
        } // for each level
    } // omp parallel
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void forwardSolve(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule,
    const int *perm)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        const int ntasks = schedule.ntasks;
        const short *nparents = schedule.nparentsForward;
        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        int nPerThread = (ntasks + nthreads - 1) / nthreads;
        int nBegin = min(nPerThread * tid, ntasks);
        int nEnd = min(nBegin + nPerThread, ntasks);

        volatile int *taskFinished = schedule.taskFinished; // 使用volatile 关键字确保每次都是从内存中读取最新的值
        int **parents = schedule.parentsForward;

        memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

        synk::Barrier::getInstance()->wait(tid);

        for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task)
        {
            SPMP_LEVEL_SCHEDULE_WAIT; // 这里是按照nparent 都求解出来才能继续这个task的求解，使用的是intrinsic中的mm_pause自旋锁，提升了性能

            for (int i = taskBoundaries[task]; i < taskBoundaries[task + 1]; ++i)
            {
                int row = perm[i] + base;
                double sum = b[row];
                for (int j = rowptr[row]; j < rowptr[row + 1]; ++j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[row] = sum * idiag[row];
            }

            SPMP_LEVEL_SCHEDULE_NOTIFY; // 将对应的task标记为 1，这里的p2p类似于sync-free
        } // for each task
    } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void backwardSolve(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule,
    const int *perm)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        const int ntasks = schedule.ntasks;
        const short *nparents = schedule.nparentsBackward;
        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        int nPerThread = (ntasks + nthreads - 1) / nthreads;
        int nBegin = min(nPerThread * tid, ntasks);
        int nEnd = min(nBegin + nPerThread, ntasks);

        volatile int *taskFinished = schedule.taskFinished;
        int **parents = schedule.parentsBackward;

        memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

        synk::Barrier::getInstance()->wait(tid);

        for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task)
        {
            SPMP_LEVEL_SCHEDULE_WAIT;

            for (int i = taskBoundaries[task + 1] - 1; i >= taskBoundaries[task]; --i)
            {
                int row = perm[i] + base;
                double sum = b[row];
                for (int j = rowptr[row + 1] - 1; j >= rowptr[row]; --j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[row] = sum;
            }

            SPMP_LEVEL_SCHEDULE_NOTIFY;
        } // for each task
    } // omp parallel
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization. Matrix is reordered.
 */
void forwardSolveWithBarrierAndReorderedMatrix(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task)
        {
            for (int i = taskBoundaries[task] + base; i < taskBoundaries[task + 1] + base; ++i)
            {
                double sum = b[i];
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[i] = sum * idiag[i];
            } // for each row
            synk::Barrier::getInstance()->wait(tid);
        } // for each level
    } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and barrier synchronization. Matrix is reordered.
 */
void backwardSolveWithBarrierAndReorderedMatrix(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();

        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task)
        {
            for (int i = taskBoundaries[task + 1] - 1 + base; i >= taskBoundaries[task] + base; --i)
            {
                double sum = b[i];
                for (int j = rowptr[i + 1] - 1; j >= rowptr[i]; --j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[i] = sum;
            } // for each row
            synk::Barrier::getInstance()->wait(tid);
        } // for each level
    } // omp parallel
}

/**
 * Forward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void forwardSolveWithReorderedMatrix(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        const int ntasks = schedule.ntasks;
        const short *nparents = schedule.nparentsForward;
        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        int nPerThread = (ntasks + nthreads - 1) / nthreads;
        int nBegin = min(nPerThread * tid, ntasks);
        int nEnd = min(nBegin + nPerThread, ntasks);

        volatile int *taskFinished = schedule.taskFinished;
        int **parents = schedule.parentsForward;

        memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

        synk::Barrier::getInstance()->wait(tid);

        for (int task = threadBoundaries[tid]; task < threadBoundaries[tid + 1]; ++task)
        {
            SPMP_LEVEL_SCHEDULE_WAIT;

            for (int i = taskBoundaries[task] + base; i < taskBoundaries[task + 1] + base; ++i)
            {
                double sum = b[i];
                for (int j = rowptr[i]; j < rowptr[i + 1]; ++j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[i] = sum * idiag[i];
            }

            SPMP_LEVEL_SCHEDULE_NOTIFY;
        } // for each task
    } // omp parallel
}

/**
 * Backward sparse triangular solver parallelized with level scheduling
 * and point-to-point synchronization
 */
void backwardSolveWithReorderedMatrix(
    const CSR &A, double y[], const double b[],
    const LevelSchedule &schedule)
{
    ADJUST_FOR_BASE;

#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        const int ntasks = schedule.ntasks;
        const short *nparents = schedule.nparentsBackward;
        const vector<int> &threadBoundaries = schedule.threadBoundaries;
        const vector<int> &taskBoundaries = schedule.taskBoundaries;

        int nPerThread = (ntasks + nthreads - 1) / nthreads;
        int nBegin = min(nPerThread * tid, ntasks);
        int nEnd = min(nBegin + nPerThread, ntasks);

        volatile int *taskFinished = schedule.taskFinished;
        int **parents = schedule.parentsBackward;

        memset((char *)(taskFinished + nBegin), 0, (nEnd - nBegin) * sizeof(int));

        synk::Barrier::getInstance()->wait(tid);

        for (int task = threadBoundaries[tid + 1] - 1; task >= threadBoundaries[tid]; --task)
        {
            SPMP_LEVEL_SCHEDULE_WAIT;

            for (int i = taskBoundaries[task + 1] - 1 + base; i >= taskBoundaries[task] + base; --i)
            {
                double sum = b[i];
                for (int j = rowptr[i + 1] - 1; j >= rowptr[i]; --j)
                {
                    sum -= values[j] * y[colidx[j]];
                }
                y[i] = sum;
            }

            SPMP_LEVEL_SCHEDULE_NOTIFY;
        } // for each task
    } // omp parallel
}


void write_res_to_csv(const char *filename,
                      const char *matrix_name,
                      const int nrow,
                      const int nnzA,
                      const int levels,
                      const int p2pNum,
                      const int p2pNumWithRT,
                      const double fwd_ref_gflops,
                      const double fwd_ref_gbps,
                      const double fwd_barrier_gflops,
                      const double fwd_barrier_gbps,
                      const double fwd_p2p_gflops,
                      const double fwd_p2p_gbps,
                      const double fwd_p2p_tr_red_gflops,
                      const double fwd_p2p_tr_red_gbps,
                      const double fwd_barrier_perm_gflops,
                      const double fwd_barrier_perm_gbps,
                      const double fwd_p2p_perm_gflops,
                      const double fwd_p2p_perm_gbps,
                      const double fwd_p2p_tr_red_perm_gflops,
                      const double fwd_p2p_tr_red_perm_gbps)

{
    char header[8192] = {'\0'};
    strcpy(header, "matrix");
    strcat(header, ",test_timestamp");
    strcat(header, ",nrow");
    strcat(header, ",nnzA");
    strcat(header, ",levels");
    strcat(header, ",p2pNum");
    strcat(header, ",p2pNumWithRT");
    strcat(header, ",fwd_ref_gflops");
    strcat(header, ",fwd_ref_gbps");
    strcat(header, ",fwd_barrier_gflops");
    strcat(header, ",fwd_barrier_gbps");
    strcat(header, ",fwd_p2p_gflops");
    strcat(header, ",fwd_p2p_gbps");
    strcat(header, ",fwd_p2p_tr_red_gflops");
    strcat(header, ",fwd_p2p_tr_red_gbps");
    strcat(header, ",fwd_barrier_perm_gflops");
    strcat(header, ",fwd_barrier_perm_gbps");
    strcat(header, ",fwd_p2p_perm_gflops");
    strcat(header, ",fwd_p2p_perm_gbps");
    strcat(header, ",fwd_p2p_tr_red_perm_gflops");
    strcat(header, ",fwd_p2p_tr_red_perm_gbps");
    char *content = (char *)g_malloc(8192 * sizeof(char));
    memset(content, 0, 8192 * sizeof(char));
    char cache[256] = {'\0'};
    // strcat(header, ",test_timestamp");
    char time_str[64] = {'\0'};
    time_t now_time;
    struct tm *p;
    time(&now_time);
    p = gmtime(&now_time);
    sprintf(time_str, "%d-%02d-%02d %02d:%02d:%02d", 1900 + p->tm_year, 1 + p->tm_mon, p->tm_mday, 8 + p->tm_hour, p->tm_min, p->tm_sec);
    sprintf(content, "%s,%s,%ld,%ld,%ld,%ld,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
            matrix_name, time_str, nrow, nnzA, levels,p2pNum, p2pNumWithRT,
            fwd_ref_gflops, fwd_ref_gbps,
            fwd_barrier_gflops, fwd_barrier_gbps,
            fwd_p2p_gflops, fwd_p2p_gbps,
            fwd_p2p_tr_red_gflops, fwd_p2p_tr_red_gbps,
            fwd_barrier_perm_gflops, fwd_barrier_perm_gbps,
            fwd_p2p_perm_gflops, fwd_p2p_perm_gbps,
            fwd_p2p_tr_red_perm_gflops, fwd_p2p_tr_red_perm_gbps);

    
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
// #pragma omp parallel for
//     for (int i = 0; i < count; i++)
//     {
//         coo->col_idx[i]++;
//         coo->row_idx[i]++;
//     }

    free(rowidx);
    free(colidx);
    free(values);
}

void SpMP_kernel(const char *matrix_name, const char *csv_filename, int num_test)
{
    // printf("file path:%s\n", matrix_name);
    int m, n;
    int nnz;
    int isSymmetric;
    int *cooRowIdx;
    int *cooColIdx;
    double *cooVal;

    int res = mmio_allinone_coo(&m, &n, &nnz, &isSymmetric,
                            &cooRowIdx, &cooColIdx, &cooVal,
                            matrix_name);
    printf("%s: %s\n", matrix_name, res == 0 ? "read mtx success!" : "read mtx error!");
    COO_MTX cooA;
    cooA.row_idx = cooRowIdx;
    cooA.col_idx = cooColIdx;
    cooA.val = cooVal;
    cooA.ncol = n;
    cooA.nrow = m;
    cooA.nnz = nnz;
    printf("original row, col, nnz: %d %d %d\n", m, n, nnz);

    pre_convert_to_lower_triangular(&cooA);
    make_full(&cooA);       // zero-based or one-based is defined by implementation
    m = cooA.nrow;
    nnz = cooA.nnz;

    CSR *A = new CSR(m, n, cooA.nnz);
    COO *coo = new COO();
    coo->m = m;
    coo->n = n;
    coo->nnz = nnz;
    coo->colidx = cooA.col_idx;
    coo->rowidx = cooA.row_idx;
    coo->values = cooA.val;

    dcoo2csr(A, coo);   // zero-based by us modification
    free(cooA.row_idx);
    free(cooA.col_idx);
    free(cooA.val);

    CSR *L = new CSR, *U = new CSR;
    splitLU(*A, L, U); // 这个函数出现了问题

    int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
    // 这个函数用以测试矩阵是否是对称的，并且获得额外的信息， 是否对称，并生成其对称形式的非零模式。
    // 如果矩阵是对称的，函数将释放所有已分配的内存并返回 true。
    // 如果矩阵不是对称的，函数将继续构造完整的对称非零模式，并返回 false
    // 构造完整对称非零模式，是将非对称矩阵构造为对称矩阵的模式，如果(i, j)的对称位置(j, i)不存在非零元就构造出这个非零元，如果存在不进行处理。
    bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);
    printf("issymmetric: %s\n", isSymmetric ? "true" : "false");

    // Todo: 下面的三行出现了segmentation error
    // int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
    // bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);
    // printf("issymmetric: %s\n", isSymmetric ? "true" : "false");    // 这里输出了false
    // fflush(stdout);


    // Todo: 下面的三行出现了segmentation error

    LevelSchedule *barrierSchedule = new LevelSchedule;
    barrierSchedule->useBarrier = true;
    barrierSchedule->transitiveReduction = false;
    if (wasSymmetric)
    {
        FREE(symRowPtr);
        FREE(symColIdx);
        FREE(symDiagPtr);
        FREE(symExtPtr);

        barrierSchedule->constructTaskGraph(*A);
    }
    else
    {
        barrierSchedule->constructTaskGraph(
            A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
            PrefixSumCostFunction(symRowPtr));
    }

    LevelSchedule *p2pSchedule = new LevelSchedule;
    p2pSchedule->transitiveReduction = false;
    if (wasSymmetric)
    {
        p2pSchedule->constructTaskGraph(*A);
    }
    else
    {
        p2pSchedule->constructTaskGraph(
            A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
            PrefixSumCostFunction(symRowPtr));
    }

    LevelSchedule *p2pScheduleWithTransitiveReduction = new LevelSchedule;
    if (wasSymmetric)
    {
        p2pScheduleWithTransitiveReduction->constructTaskGraph(*A);
    }
    else
    {
        p2pScheduleWithTransitiveReduction->constructTaskGraph(
            A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
            PrefixSumCostFunction(symRowPtr));

        FREE(symRowPtr);
        FREE(symColIdx);
        FREE(symDiagPtr);
        FREE(symExtPtr);
    }

    printf("parallelism %f\n", (double)A->m / (barrierSchedule->levIndices.size() - 1));
    int nlevels = barrierSchedule->levIndices.size() - 1;
    assert(barrierSchedule->levIndices.size() == p2pSchedule->levIndices.size());
    assert(barrierSchedule->levIndices.size() == p2pScheduleWithTransitiveReduction->levIndices.size());

    int p2pTaskNum = p2pSchedule->ntasks;
    int p2pTaskNumTR = p2pScheduleWithTransitiveReduction->ntasks;

    int p2pNum = 0, p2pNumWithRT = 0;
    #pragma omp parallel for reduction(+ : p2pNum)
    for(int t=0; t < p2pTaskNum; t++)
    {
        p2pNum += p2pSchedule->nparentsForward[t];
    }

    #pragma omp parallel for reduction(+ : p2pNumWithRT)
    for(int t=0; t < p2pTaskNumTR; t++)
    {
        p2pNumWithRT += p2pScheduleWithTransitiveReduction->nparentsForward[t];
    }

    /////////////////////////////////////////////////////////////////////////////
    // Reorder matrix
    /////////////////////////////////////////////////////////////////////////////

    const int *perm = p2pScheduleWithTransitiveReduction->origToThreadContPerm;
    const int *invPerm = p2pScheduleWithTransitiveReduction->threadContToOrigPerm;
    assert(isPerm(perm, A->m));
    assert(isPerm(invPerm, A->m));

    CSR *LPerm = A->permute(perm, invPerm, true);
    // CSR *UPerm = U->permute(perm, invPerm, true);

    /////////////////////////////////////////////////////////////////////////////
    // Allocate vectors
    /////////////////////////////////////////////////////////////////////////////

    double *b = MALLOC(double, A->m);
#pragma omp parallel for
    for (int i = 0; i < A->m; i++)
        b[i] = i;

    double *y = MALLOC(double, A->m);
    double *x = MALLOC(double, A->m);

    double flop, byte;
    for (int i = 0; i < 2; ++i)
    {
        int nnz = L->rowptr[A->m] + A->m;
        flop = 2 * nnz;
        byte = (nnz) * (sizeof(double) + sizeof(int)) + A->m * sizeof(int);
    }

    // allocate a large buffer to flush out cache
    bufToFlushLlc = (double *)_mm_malloc(LLC_CAPACITY, 64);

    // ///////////////////////////////////////////////////////////////////////////
    // sparse triangular solver w/o reordering
    // ///////////////////////////////////////////////////////////////////////////

    int REPEAT = num_test;
    double timesForward[REPEAT], timesBackward[REPEAT];
    map<string, double> resMap;
    for (int o = REFERENCE; o <= P2P_WITH_TRANSITIVE_REDUCTION; ++o)
    {
        SynchronizationOption option = (SynchronizationOption)o;

        for (int i = 0; i < REPEAT; ++i)
        {
            flushLlc();

            initializeX(x, A->m);
            initializeX(y, A->m);

            // double t = omp_get_wtime();
            double t = getCurrentTimeMilli();

            switch (option)
            {
            case REFERENCE:
                forwardSolveRef(*L, y, b);
                break;
            case BARRIER:
                forwardSolveWithBarrier(*L, y, b, *barrierSchedule, invPerm);
                break;
            case P2P:
                forwardSolve(*L, y, b, *p2pSchedule, invPerm);
                break;
            case P2P_WITH_TRANSITIVE_REDUCTION:
                forwardSolve(*L, y, b, *p2pScheduleWithTransitiveReduction, invPerm);
                break;
            default:
                assert(false);
                break;
            }

            timesForward[i] = getCurrentTimeMilli() - t;
            if (i == REPEAT - 1)
            {
                string res_key;
                for (int j = 0; j < 1; ++j)
                {
                    printf(0 == j ? "fwd_" : "bwd_");
                    switch (option)
                    {
                    case REFERENCE:
                        printf("ref\t\t\t");
                        res_key = "fwd_ref";
                        break;
                    case BARRIER:
                        printf("barrier\t\t");
                        res_key = "fwd_barrier";
                        break;
                    case P2P:
                        printf("p2p\t\t\t");
                        res_key = "fwd_p2p";
                        break;
                    case P2P_WITH_TRANSITIVE_REDUCTION:
                        printf("p2p_tr_red\t\t");
                        res_key = "fwd_p2p_tr_red";
                        break;
                    default:
                        assert(false);
                        break;
                    }
                    printEfficiency(
                        0 == j ? timesForward : timesBackward, REPEAT, flop, byte);
                }

                sort(timesForward, timesForward + REPEAT);
   
                double time_cost = timesForward[REPEAT / 2];
                resMap.insert(make_pair(res_key + "_gflops", flop / time_cost / 1e6));
                resMap.insert(make_pair(res_key + "_gbps", byte / time_cost / 1e6));
                correctnessCheck(A, x);
            }
        } // for each iteration
    } // for each option

    // ///////////////////////////////////////////////////////////////////////////
    // // sparse triangular solver w/ reordering
    // ///////////////////////////////////////////////////////////////////////////

    double *bPerm = getReorderVector(b, perm, A->m);
    double *tempVector = MALLOC(double, A->m);

    for (int o = BARRIER; o <= P2P_WITH_TRANSITIVE_REDUCTION; ++o)
    {
        SynchronizationOption option = (SynchronizationOption)o;

        for (int i = 0; i < REPEAT; ++i)
        {
            flushLlc();

            initializeX(x, A->m);
            initializeX(y, A->m);
            reorderVector(x, tempVector, perm, A->m);
            reorderVector(y, tempVector, perm, A->m);

            double t = omp_get_wtime();

            switch (option)
            {
            case BARRIER:
                forwardSolveWithBarrierAndReorderedMatrix(
                    *LPerm, y, bPerm, *barrierSchedule);
                break;
            case P2P:
                forwardSolveWithReorderedMatrix(
                    *LPerm, y, bPerm, *p2pSchedule);
                break;
            case P2P_WITH_TRANSITIVE_REDUCTION:
                forwardSolveWithReorderedMatrix(
                    *LPerm, y, bPerm, *p2pScheduleWithTransitiveReduction);
                break;
            default:
                assert(false);
                break;
            }

            timesForward[i] = omp_get_wtime() - t;
            if (i == REPEAT - 1)
            {
                string res_key;
                for (int j = 0; j < 1; ++j)
                {
                    printf(0 == j ? "fwd_" : "bwd_");
                    switch (option)
                    {
                    case BARRIER:
                        printf("barrier_perm\t");
                        res_key = "fwd_barrier_perm";
                        break;
                    case P2P:
                        printf("p2p_perm\t\t");
                        res_key = "fwd_p2p_perm";
                        break;
                    case P2P_WITH_TRANSITIVE_REDUCTION:
                        printf("p2p_tr_red_perm\t");
                        res_key = "fwd_p2p_tr_red_perm";
                        break;
                    default:
                        assert(false);
                        break;
                    }
                    printEfficiency(
                        0 == j ? timesForward : timesBackward, REPEAT, flop, byte);
                }
                sort(timesForward, timesForward + REPEAT);

                double time_cost = timesForward[REPEAT / 2];
                resMap.insert(make_pair(res_key + "_gflops", flop / time_cost / 1e9));
                resMap.insert(make_pair(res_key + "_gbps", byte / time_cost / 1e9));
                /**
                 * 这里输出的CSV文件中
                 */
                reorderVector(x, tempVector, invPerm, A->m);
                correctnessCheck(A, x);
            }
        }
    }

    /**
     * 这里将总的结果输出到CSV文件中
     */
    write_res_to_csv(csv_filename, matrix_name, L->m, L->getNnz() + L->m, nlevels,p2pNum, p2pNumWithRT,
                     resMap["fwd_ref_gflops"],
                     resMap["fwd_ref_gbps"],
                     resMap["fwd_barrier_gflops"],
                     resMap["fwd_barrier_gbps"],
                     resMap["fwd_p2p_gflops"],
                     resMap["fwd_p2p_gbps"],
                     resMap["fwd_p2p_tr_red_gflops"],
                     resMap["fwd_p2p_tr_red_gbps"],
                     resMap["fwd_barrier_perm_gflops"],
                     resMap["fwd_barrier_perm_gbps"],
                     resMap["fwd_p2p_perm_gflops"],
                     resMap["fwd_p2p_perm_gbps"],
                     resMap["fwd_p2p_tr_red_perm_gflops"],
                     resMap["fwd_p2p_tr_red_perm_gbps"]);

#ifdef MKL

    /////////////////////////////////////////////////////////////////////////////
    // Inspector-Executor interface in MKL 11.3+
    /////////////////////////////////////////////////////////////////////////////

    sparse_matrix_t mklA;
    sparse_status_t stat = mkl_sparse_d_create_csr(
        &mklA,
        SPARSE_INDEX_BASE_ZERO, A->m, A->n,
        A->rowptr, A->rowptr + 1,
        A->colidx, A->values);

    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to create mkl csr\n");
        return -1;
    }

    matrix_descr descL;
    descL.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descL.mode = SPARSE_FILL_MODE_LOWER;
    descL.diag = SPARSE_DIAG_NON_UNIT;

    sparse_matrix_t mklL, mklU;
    stat = mkl_sparse_copy(mklA, descL, &mklL);
    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to create mkl csr lower\n");
        return -1;
    }

    stat = mkl_sparse_set_sv_hint(
        mklL, SPARSE_OPERATION_NON_TRANSPOSE, descL, REPEAT);

    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to set sv hint\n");
        return -1;
    }

    stat = mkl_sparse_optimize(mklL);

    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to sparse optimize\n");
        return -1;
    }

    matrix_descr descU;
    descU.type = SPARSE_MATRIX_TYPE_TRIANGULAR;
    descU.mode = SPARSE_FILL_MODE_UPPER;
    descU.diag = SPARSE_DIAG_UNIT;

    stat = mkl_sparse_copy(mklA, descU, &mklU);
    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to create mkl csr upper\n");
        return -1;
    }

    stat = mkl_sparse_set_sv_hint(
        mklU, SPARSE_OPERATION_NON_TRANSPOSE, descU, REPEAT);

    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to set sv hint\n");
        return -1;
    }

    stat = mkl_sparse_optimize(mklU);

    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to sparse optimize\n");
        return -1;
    }

    for (int i = 0; i < REPEAT; ++i)
    {
        flushLlc();

        initializeX(x, A->m);
        initializeX(y, A->m);

        double t = omp_get_wtime();

        mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklL, descL, b, y);

        timesForward[i] = omp_get_wtime() - t;
        t = omp_get_wtime();

        mkl_sparse_d_trsv(SPARSE_OPERATION_NON_TRANSPOSE, 1, mklU, descU, y, x);

        timesBackward[i] = omp_get_wtime() - t;

        if (i == REPEAT - 1)
        {
            for (int j = 0; j < 2; ++j)
            {
                printf(0 == j ? "fwd_" : "bwd_");
                printf("mkl\t\t\t");
                printEfficiency(
                    0 == j ? timesForward : timesBackward, REPEAT, flop, byte);
            }

            correctnessCheck(A, x);
        }
    }

    stat = mkl_sparse_destroy(mklA);
    if (SPARSE_STATUS_SUCCESS != stat)
    {
        fprintf(stderr, "Failed to destroy mkl csr\n");
        return -1;
    }

#endif

    delete barrierSchedule;
    delete p2pSchedule;
    delete p2pScheduleWithTransitiveReduction;

    delete A;
    delete L;
    delete U;

    delete LPerm;
    // delete UPerm;

    FREE(b);
    FREE(y);
    FREE(bPerm);
    FREE(tempVector);

    // delete L;
    // delete U;


    synk::Barrier::deleteInstance();

//     return 0;
fflush(stdout);
}

void run(const char *matrix_name,
        //  const char *matrix_name,
         const char *csv_res_file, int num_test)
{
    printf("###################### %s #######################\n", matrix_name);
    SpMP_kernel(matrix_name, csv_res_file, num_test);
    // SpMP_kernel();
}

int main(int argc, char *argv[])
{
    int nthreads = atoi(argv[3]);
    int num_test = atoi(argv[4]);
     printf("nthread: %d\n",nthreads);
    omp_set_num_threads(nthreads);

    char matrix_file[1024] = {'\0'};
    for (int i = 0; i < argc; i++)
    {
        printf("arg[%d]: %s\n", i, argv[i]);
    }
    // int write_res_choice = atoi(argv[3]);
    strcpy(matrix_file, argv[1]);
    // printf("matrix list file: %s\n", matrix_file);
    char csv_res_file[1024] = "\0";
    strcpy(csv_res_file, argv[2]);

    // char matrix_source_dir[1024] = "/fs2/home/nudt_liujie/hch/matrix_mtx/";
    // printf("read matrix list end!!\n");

    // char *matrix_s_dir_temp = (char *)g_calloc(1024, sizeof(char));
    // strcpy(matrix_s_dir_temp, matrix_source_dir);
    run(matrix_file , csv_res_file, num_test);
    // memset(matrix_s_dir_temp, 0, 1024 * sizeof(char));
    // free(matrix_s_dir_temp);
}