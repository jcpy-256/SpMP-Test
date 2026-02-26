#ifndef _UTILS_H
#define _UTILS_H
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>

#include "matrix.h"
#include "common_macros.h"


/*----------------------- time saved in all bench repeat ---------------------------*/
// 时间是dsp的 clock 数
typedef struct stage_time_cost_clk
{
    // uint64_t time_read_file;          // 读取文件
    // uint64_t time_analysis;    // 预处理
    // uint64_t time_pre_triangular; // 处理为下三角的时间
    // uint64_t time_coo2csr; //将 coo转换为csr的时间
    // uint64_t time_solve;    // 求解的时间

    // 这部分时间，每个轮次中取最大值求累加
    uint64_t time_aggreate_solution; // 聚合解的时间
   // uint64_t time_transfer_data;     // 加载数据
   // uint64_t time_calculate;         // 计算时间

   // 传输非对角快的时间
    uint64_t time_transfer_nondiag_val;
    // 搬运x的时间
    uint64_t time_transfer_x;
    // 搬运b的时间
    uint64_t time_transfer_b;
    // 传输对角块的时间
    uint64_t time_transfer_diag_data;
    // 部分和的计算时间
    uint64_t time_left_calculate;
    // 对角块求解时间
    uint64_t time_diag_calculate;
    uint64_t time_barrier;  // 同步时间
    uint64_t time_interior; // dsp端的执行clock

} STC_CLK;

/* -----------------typedef matrix struct name -------------- */
typedef struct CSR_Matrix CSR_MTX;
typedef struct COO_Matrix COO_MTX;
typedef struct ELL_Matrix ELL_MTX;
typedef struct CSC_Matrix CSC_MTX;
typedef struct BSR_Matrix BSR_MTX;


void CSR_Free(CSR_MTX csrA)
{
    FREE_PTR(csrA.col_idx);
    FREE_PTR(csrA.row_ptr);
    FREE_PTR(csrA.val);
}

/* ------------------matrix operation------------------------ */
typedef struct level_set
{
    int *level_list;             // every node level_list,  level == 1- based
    int *level_ptr;              // node number ptr in one level, length = nlevels + 1
    int *nodes_grouped_by_level; // nodes grouped by level, length = row_nums 0-based
    int nlevels;                 // level number
    int max_in;
} LEV_SET;

typedef union matrix_store_format
{
    COO_MTX *coo;
    CSR_MTX *csr;
    ELL_MTX *ELL;
    BSR_MTX *bsr;

} MTX_STORE;

typedef enum matrix_format_enum
{
    COO_ENUM,
    CSR_ENUM,
    CSC_ENUM,
    DIA_ENUM,
    ELL_ENUM,
    BSR_ENUM
} MTX_FMT_ENUM;

typedef struct matrix_info
{
    int nrow; // 存储矩阵中真实的nrow 而非cbsr或者其他方式处理过的nrow
    int ncol;
    long nnz;
} MTX_INFO;

/**
 * 可以记录使用线程的信息的传递的参数的信息。
 */
struct thsparseContext
{
    int clusterId;
    int thread_count;
    char datFilename[1024];
    char kernel[256];
    uint64_t time_cost;
    STC_CLK stc_clk;
    int rwlock_id;
    int barrier_id;
};

/* mat*/
struct thsparseSpMatDescr
{
    MTX_INFO info;
    MTX_FMT_ENUM fmt; /*matrix format*/
    MTX_STORE matrix; /* matrix store union in store_format */
};

/* vector */
struct thsparseDnVecDescr
{
    double *values;
    int len;
};
/* sptrsv description */
struct thsparseSpTrSVDescr
{
    /* data */
    LEV_SET *level_set;
    double *diagonal_blk_val;  // 存储非对角快
};

/* matrix operation */
typedef enum thsparseOperation_t
{
    NULL_OPERATION,
    TRANSPOSITION
} thsparseOperation_t;

/* algorithm */
typedef enum thsparseSpTRSVAlg_t
{
    LEVEL_SET_ALG,
    SYNC_FREE_ALG
} thsparseSpTRSVAlg_t;

typedef enum thsparseStatus_t
{
    SUCCESS,
    ERROR
} thsparseStatus_t;

typedef enum thsparseDataType_t
{
    SPARK_32F,
    SPARK_64F
} thsparseDataType_t;

// ############################ cuSPARSE API Opaque Data Structure ########################

// execition context
typedef struct thsparseContext *thsparseHandle_t;

typedef struct thsparseSpMatDescr *thsparseSpMatDescr_t;

typedef struct thsparseDnVecDescr *thsparseDnVecDescr_t;

// analsysi phrase information
typedef struct thsparseSpTrSVDescr *thsparseSpTrSVDescr_t;

// /**
//  * @brief using malloc function to allocate memory
//  * @param size: byte size needed to allocate
//  *
//  */
// inline void *g_malloc(size_t size)
// {
//     void *p = malloc(size);
//     if (p == NULL)
//     {
//         perror("memory allocation (malloc) failure!\n");
//     }
//     return (p);
// }

// /**
//  * @brief using calloc function to allocate memory and the element is initialized with 0
//  */
// inline void *g_calloc(size_t num, size_t size)
// {
//     void *p = calloc(num, size);
//     if (p == NULL)
//     {
//         perror("memory allocation (calloc) failure!\n");
//     }
//     return (p);
// }

// /**
//  * @brief After successfully applying size bytes memory, Copy size bytes data from source to dest and free source memory
//  */
// inline void g_memcpy(void **dest, void *source, size_t size)
// {
//     void *p = g_malloc(size);
//     memcpy(p, source, size);
//     *dest = p;
//     free(source);
// }

// void mtx_enum_to_str(MTX_FMT_ENUM mfe,char* type)
// {
//     // char type[8];
//     switch (mfe)
//     {
//     case COO:
//         strcpy(type, COO_FROMAT);
//         break;
//     case CSC:
//         strcpy(type, CSC_FROMAT);
//         break;
//     case CSR:
//         strcpy(type, CSR_FORMAT);
//         break;
//     case ELL:
//         strcpy(type, DIA_FROMAT);
//         break;
//     case DIA:
//         strcpy(type, DIA_FROMAT);
//         break;
//     default:
//         break;
//     }
//     printf("%s\n",type);
//     return type;
// }

// inline void prefix_sum(int *input, int length)
// {
//     if (length == 0 || length == 1)
//         return;

//     // int old_val, new_val;

//     // old_val = input[0];
//     // input[0] = 0;
//     // for (int i = 1; i < length; i++)
//     // {
//     //     new_val = input[i];
//     //     input[i] = old_val + input[i - 1];
//     //     old_val = new_val;
//     // }

//     for (int i = 1; i < length; i++)
//     {
//         input[i] = input[i - 1] + input[i];
//     }
// }

// inline uint64_t getCurrentTimeMilli()
// {
//     struct timeval time;

//     gettimeofday(&time, NULL);
//     return (uint64_t)(time.tv_sec * INT64_C(1000) + time.tv_usec / 1000);
// }

// inline int fileIsExist(const char *filename)
// {
//     return access(filename, F_OK);
// }

// inline uint64_t doubleToRawBits(double d)
// {
//     union
//     {
//         uint64_t i;
//         double f;
//     } word;
//     word.f = d;
//     return word.i;
// }

#endif