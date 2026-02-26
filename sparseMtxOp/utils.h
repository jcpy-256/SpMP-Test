/**
 * 这个文件存放共有的函数，主要是功能函数，例如获取当前的时间等。
 */

#ifndef UTILS_H
#define UTILS_H
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>


/**
 * @brief using malloc function to allocate memory
 * @param size: byte size needed to allocate
 *
 */
void *g_malloc(size_t size)
{
    void *p = malloc(size);
    if (p == NULL)
    {
        perror("memory allocation (malloc) failure!\n");
    }
    return (p);
}

/**
 * @brief using calloc function to allocate memory and the element is initialized with 0
 */
void *g_calloc(size_t num, size_t size)
{
    void *p = calloc(num, size);
    if (p == NULL)
    {
        perror("memory allocation (calloc) failure!\n");
    }
    return (p);
}

/**
 * @brief After successfully applying size bytes memory, Copy size bytes data from source to dest and free source memory
 */
void g_memcpy(void **dest, void *source, size_t size)
{
    void *p = g_malloc(size);
    memcpy(p, source, size);
    *dest = p;
    free(source);
}

/**
 * get current time in
 */

/**
 * get current time in microsecond
 */
uint64_t getCurrentTimeMicro()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return (uint64_t)(time.tv_sec * INT64_C(1000000) + time.tv_usec);
}

double getCurrentTimeMilli()
{
    struct timeval time;

    gettimeofday(&time, NULL);
    return (double)(time.tv_sec * INT64_C(1000) + 1.0 * time.tv_usec / 1000);
}
/**
 * res == 0 file exists
 */
int fileIsExist(const char *filename)
{
    return access(filename, F_OK);
}

/**
 * return 1 pass
 * retunr 
 */
int check_solutions_double(double *check, double *correct, int length, const double error)
{
    printf("#################check value!\n");
    for (int i = 0; i < length; i++)
    {
        if (fabs(check[i] - correct[i]) > error)
        {

            printf("check value:%lf, correct value: %lf\n", check[i], correct[i]);
            return -1 * i;
        }
    }
    return 1;
}

void array_init_uint64(uint64_t *array, uint64_t value, int len)
{
    for (int i = 0; i < len; i++)
    {
        array[i] = value;
    }
}

void array_init_int(int *array, int value, int len)
{
    for (int i = 0; i < len; i++)
    {
        array[i] = value;
    }
}

void array_init_double(double *array, double value, int len)
{
    for (int i = 0; i < len; i++)
    {
        array[i] = value;
    }
}

// double getTimeMilliByCLK(uint64_t clk)
// {
//     return 1.0 * clk * MILLISECOND / DSP_CLK_FREQ;
// }

void write_soulutions(const char *filename, double *solutions, const int len)
{
    FILE *fp = fopen(filename, "w+");
    if (fp == NULL)
    {
        perror(" file open error!\n");
    }
    for (int i = 0; i < len; i++)
    {
        // printf("x[%d]: %lf\n", i, solutions_blk[i]);
        fprintf(fp, "x[%d]: %lf\n", i, solutions[i]);
        // if(i == 27696){
        //     printf("27696 : %lf\n", solutions_blk[i]);
        // }
    }
    fclose(fp);
}


/**
 * 如果文件不存在则创建文件，并将文件头的内容写入，
 * 文件存在就将item写入，
 * 函数自动加入换行，
 * 
 */
void wirte_item_to_csv(const char* fileheader, const char *item, const char *filename)
{
    // printf("%s", fileheader);
    // printf("%s", item);
    // printf("file exist: %d\n", fileIsExist(filename));
    // printf("header len %d\n", strlen(fileheader));
    // printf("item len %d\n", strlen(item));
    if(fileIsExist(filename) == 0 && strlen(item) > 0) // 将文件打开后，追加一行数据
    {
        FILE *fp = fopen(filename, "a");
         if(fp == NULL)
        {
            printf("file %s open error!\n",filename);
            exit(2);
        }else{
            fprintf(fp, "%s\n",item);
        }
        fclose(fp);
    }else if(strlen(fileheader) > 0){ // 文件不存在创建文件，并将文件头写入
        FILE *fp = fopen(filename, "w+");
        if(fp == NULL)
        {
            printf("file %s open error!\n",filename);
            exit(2);
        }else{
            printf("write header and content\n");
            fprintf(fp, "%s\n",fileheader);
            fprintf(fp, "%s\n",item);
        }
        fclose(fp);
    }else{
        printf("write item to csv function parameter is error!\n");
        exit(2);
    }
}

#endif