/**
 * 2025-09-08 by yansen
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

inline double sparse_dot_product(int l1, int u1, int l2, int u2, const int *indices, const double *data)
{
    double result = 0.0;
    while (l1 < u1 && l2 < u2)
    {
        if (indices[l1] == indices[l2]) // matching column?
        {
            assert(std::isnan(data[l1]) == false);
            assert(std::isnan(data[l2]) == false);
            assert(std::isinf(data[l1]) == false);
            assert(std::isinf(data[l2]) == false);
            result += data[l1++] * data[l2++];
        }
        else if (indices[l1] < indices[l2]) // else proceed until we find matching columns
        {
            l1++;
        }
        else
        {
            l2++;
        }
    }
    return result;
}

//=============================== Up looking Looking ==============================
void spic0_csr_uL_serial(const CSR *A, double *lu)
{
    const int n = A->n;
    const int *rowptr = A->rowptr;
    const int *colidx = A->colidx;
    const int *diagptr = A->diagptr;
    const double *values = A->values;

    // bool finish_flag = true;
    for (int i = 0; i < n; ++i)
    {
        for (int k = rowptr[i]; k < diagptr[i]; ++k)
        {
            const int j = colidx[k]; // column
            double dp = sparse_dot_product(
                rowptr[i], diagptr[i], // i-th row minus diagonal
                rowptr[j], diagptr[j], // j-th row minus diagonal
                colidx, lu);

            const double A_ij = values[k];

            // below diagonal?
            const double L_jj = lu[diagptr[j]]; // diagonal is last entry of j-th row
            // assert(fabs(L_jj) > 1.0 - 16);
            //    /
            assert(L_jj > 1.0e-10);
            assert(std::isnan(dp) == false);
            assert(std::isnan(A_ij) == false);
            assert(std::isnan(A_ij - dp) == false);
            assert(std::isnan(L_jj) == false);

            assert(std::isinf(dp) == false);
            assert(std::isinf(A_ij) == false);
            assert(std::isinf(A_ij - dp) == false);
            assert(std::isinf(L_jj) == false);
            lu[k] = (A_ij - dp) / L_jj;
            assert(std::isnan(lu[k]) == false);
            assert(std::isinf(lu[k]) == false);
            // }/
        }

        // update diagonal element
        double dv = sparse_dot_product(
            rowptr[i], diagptr[i], // i-th row minus diagonal
            rowptr[i], diagptr[i], // j-th row minus diagonal
            colidx, lu);
        const double A_ii = values[diagptr[i]];
        // assert(A_ii - dv >= 0);
        if (A_ii - dv > 0)
        {
            assert(std::isnan(A_ii - dv) == false);
            assert(std::isinf(A_ii - dv) == false);
            lu[diagptr[i]] = std::sqrt(A_ii - dv);
            assert(std::isnan(lu[diagptr[i]]) == false);
            assert(std::isinf(lu[diagptr[i]]) == false);
        }
    }
}

// parallel ilu0
void ic0(const CSR &A, double *lu, const LevelSchedule &schedule)
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

            for (int t = taskBoundaries[task]; t < taskBoundaries[task + 1]; ++t)
            {
                int i = perm[t] + base;

                for (int k = rowptr[i]; k < diagptr[i]; k++)
                {
                    const int j = colidx[k]; // column
                    double dp = sparse_dot_product(
                        rowptr[i], diagptr[i], // i-th row minus diagonal
                        rowptr[j], diagptr[j], // j-th row minus diagonal
                        colidx, lu);
                    const double A_ij = lu[k];
                    // below diagonal
                    const double L_jj = lu[diagptr[j]]; // diagonal element
                    assert(L_jj > 1.0e-10);
                    assert(std::isnan(dp) == false);
                    assert(std::isnan(A_ij) == false);
                    assert(std::isnan(A_ij - dp) == false);
                    assert(std::isnan(L_jj) == false);

                    assert(std::isinf(dp) == false);
                    assert(std::isinf(A_ij) == false);
                    assert(std::isinf(A_ij - dp) == false);
                    assert(std::isinf(L_jj) == false);
                    lu[k] = (A_ij - dp) / L_jj;
                    assert(std::isnan(lu[k]) == false);
                    assert(std::isinf(lu[k]) == false);
                } // for: rowptr[i] to diagptr[i]

                // update diagonal element
                double dv = sparse_dot_product(
                    rowptr[i], diagptr[i],
                    rowptr[i], diagptr[i],
                    colidx, lu);
                const double A_ii = lu[diagptr[i]];
                if (A_ii - dv > 0)
                {
                    assert(std::isnan(A_ii - dv) == false);
                    assert(std::isinf(A_ii - dv) == false);
                    lu[diagptr[i]] = std::sqrt(A_ii - dv);
                    assert(std::isnan(lu[diagptr[i]]) == false);
                    assert(std::isinf(lu[diagptr[i]]) == false);
                }
            } // for each row

            SPMP_LEVEL_SCHEDULE_NOTIFY;
        } // for each level
    } // omp parallel
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

void IC0_kernel(const char *matrix_name, const char *csv_filename, int num_test)
{
    // printf("file path:%s\n", matrix_name);
    int m, n;
    int nnz;
    CSR *A = new CSR(matrix_name);
    A->make0BasedIndexing();
    m = A->m, n = A->n;
    nnz = A->getNnz();

    // int num_test = 100;

    // ************************** Serial IC0 Test ******************************
    std::vector<double> timeVec(num_test, 0.0);
    double *lu_serial = MALLOC(double, nnz);

    double a = getCurrentTimeMilli();
    for (int i = 0; i < num_test; i++)
    {
        // copy value to lu_serial
        copyVector(lu_serial, A->values, nnz);
        double t1 = getCurrentTimeMilli();
        spic0_csr_uL_serial(A, lu_serial);
        timeVec[i] = getCurrentTimeMilli() - t1;
    }

    printf("serial time using is %lf\n", getCurrentTimeMilli() - a);
    fflush(stdout);
    // get median
    std::sort(timeVec.begin(), timeVec.end());
    double serial_time = num_test % 2 == 1 ? timeVec[(num_test) / 2] : (timeVec[(num_test + 1) / 2] + timeVec[(num_test) / 2]) / 2;

    // ************************** p2p IC0 Test ******************************
    LevelSchedule p2pScheduleWithTransitiveReduction;
    int *symRowPtr = NULL, *symColIdx = NULL, *symDiagPtr = NULL, *symExtPtr = NULL;
    bool wasSymmetric = getSymmetricNnzPattern(A, &symRowPtr, &symDiagPtr, &symExtPtr, &symColIdx);
    // printf("issymmetric: %s\n", isSymmetric ? "true" : "false");
    // p2pScheduleWithTransitiveReduction.constructTaskGraph(*A);
    if (wasSymmetric)
    {
        p2pScheduleWithTransitiveReduction.constructTaskGraph(*A);
    }
    else
    {
        p2pScheduleWithTransitiveReduction.constructTaskGraph(
            A->m, symRowPtr, symDiagPtr, symExtPtr, symColIdx,
            PrefixSumCostFunction(symRowPtr));
    }

    FREE(symRowPtr);
    FREE(symColIdx);
    FREE(symDiagPtr);
    FREE(symExtPtr);
    printf("parallelism = %g\n", (double)A->m / (p2pScheduleWithTransitiveReduction.levIndices.size() - 1));
    int nlevels = p2pScheduleWithTransitiveReduction.levIndices.size() - 1;
    double *lu_p2p = MALLOC(double, A->getNnz());
    for (int i = 0; i < num_test; i++)
    {
        // copy value to lu_serial
        copyVector(lu_p2p, A->values, nnz);
        double t1 = getCurrentTimeMilli();
        ic0(*A, lu_p2p, p2pScheduleWithTransitiveReduction);
        timeVec[i] = getCurrentTimeMilli() - t1;
        // printf("turn: %d\n", i);
    }
    // get median
    std::sort(timeVec.begin(), timeVec.end());
    double p2p_time = num_test % 2 == 1 ? timeVec[(num_test) / 2] : (timeVec[(num_test + 1) / 2] + timeVec[(num_test) / 2]) / 2;

    write_res_to_csv(csv_filename, matrix_name, num_test, A->m, A->getNnz(), nlevels, serial_time, p2p_time);
    printf("serial time: %lf, p2p time: %lf\n", serial_time, p2p_time);
    bool check = correctnessCheck(lu_p2p, lu_serial, A->getNnz(), 1.0e-10);
    if (check)
    {
        printf("IC check pass!\n");
    }
    else
    {
        printf("IC check not pass!\n");
    }
    FREE(lu_serial);
    FREE(lu_p2p);
    delete A;
}

void run(         const char *matrix_name,
         const char *csv_res_file, int num_test)
{
    // char time_str[64] = {'\0'};
    // time_t now_time;
    // struct tm *p;
    // time(&now_time);
    // p = gmtime(&now_time);
    // sprintf(time_str, "%d-%d-%d %d:%d:%d %s.txt", 1900 + p->tm_year, 1 + p->tm_mon, p->tm_mday, 8 + p->tm_hour, p->tm_min, p->tm_sec, matrix_name);
    // // strcat(performance_file, time_str);
    // char mtx_file[1024] = {'\0'};
    // strcat(matrix_source_dir, matrix_name);
    // strcpy(mtx_file, matrix_source_dir);
    // strcat(mtx_file, ".mtx");
    // printf("mtx file: %s\n", mtx_file);
    printf("###################### %s #######################\n", matrix_name);

    // printf("performance file: %s\n", performance_file);
    // sptrsv_cuda_kernel_cusparse(matrix_source_dir, file_base, res_base, matrix_name, performance_file, csv_res_file, write_res_choice, cusparse_version);
    IC0_kernel( matrix_name, csv_res_file, num_test);
}

int main(int argc, char *argv[])
{
    int nthreads = atoi(argv[3]);
    int num_test = atoi(argv[4]);
    printf("nthread: %d\n",nthreads);
    omp_set_num_threads(nthreads);

    char matrix_name[1024] = {'\0'};

    strcpy(matrix_name, argv[1]);
    // printf("matrix list file: %s\n", matrix_name);
    char csv_res_file[1024] = "\0";
    strcpy(csv_res_file, argv[2]);

    // char matrix_source_dir[1024] = "/fs2/home/nudt_liujie/hch/matrix_mtx/";
    // printf("read matrix list end!!\n");

    // char *matrix_s_dir_temp = (char *)g_calloc(1024, sizeof(char));
    // strcpy(matrix_s_dir_temp, matrix_source_dir);
    run(matrix_name, csv_res_file, num_test);
    // printf("end run!\n");
    // fflush(stdout);
    // memset(matrix_s_dir_temp, 0, 1024 * sizeof(char));
    // free(matrix_s_dir_temp);
}