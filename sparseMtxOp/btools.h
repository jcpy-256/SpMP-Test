#ifndef BTOOLS_H
#define BTOOLS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include "common.h"

/**
 * quick sort [left, right] from small to large
 */
void quicksort(int *idx, double *w, int left, int right)
{
    if (left >= right)
        return;

    std::swap(idx[left], idx[left + (right - left) / 2]);
    std::swap(w[left], w[left + (right - left) / 2]);

    int last = left;
    for (int i = left + 1; i <= right; i++)
    {
        if (idx[i] < idx[left])
        {
            ++last;
            std::swap(idx[last], idx[i]);
            std::swap(w[last], w[i]);
        }
    }

    std::swap(idx[left], idx[last]);
    std::swap(w[left], w[last]);

    quicksort(idx, w, left, last - 1);
    quicksort(idx, w, last + 1, right);
}


void exclusive_scan(int *input, int length)
{
    if (length == 0 || length == 1)
        return;

    int old_val, new_val;

    old_val = input[0];
    input[0] = 0;
    for (int i = 1; i < length; i++)
    {
        new_val = input[i];
        input[i] = old_val + input[i - 1];
        old_val = new_val;
    }
}

void swap_key(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void swap_val(MAT_VAL_TYPE *a, MAT_VAL_TYPE *b)
{
    MAT_VAL_TYPE tmp = *a;
    *a = *b;
    *b = tmp;
}

// quick sort key-value pair (child function)
int partition_key_val_pair(int *key, MAT_VAL_TYPE *val, int length, int pivot_index)
{
    int i = 0;
    int small_length = pivot_index;

    int pivot = key[pivot_index];
    swap_key(&key[pivot_index], &key[pivot_index + (length - 1)]);
    swap_val(&val[pivot_index], &val[pivot_index + (length - 1)]);

    for (; i < length; i++)
    {
        if (key[pivot_index + i] < pivot)
        {
            swap_key(&key[pivot_index + i], &key[small_length]);
            swap_val(&val[pivot_index + i], &val[small_length]);
            small_length++;
        }
    }

    swap_key(&key[pivot_index + length - 1], &key[small_length]);
    swap_val(&val[pivot_index + length - 1], &val[small_length]);

    return small_length;
}

// quick sort key-value pair (main function)
void quick_sort_key_val_pair(int *key, MAT_VAL_TYPE *val, int length)
{
    if (length == 0 || length == 1)
        return;

    int small_length = partition_key_val_pair(key, val, length, 0);
    quick_sort_key_val_pair(key, val, small_length);
    quick_sort_key_val_pair(&key[small_length + 1], &val[small_length + 1], length - small_length - 1);
}

void swap_int(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

int choose_pivot(int i, int j)
{
    return (i + j) / 2;
}

/**
 * 根据key 对value 进行快速排序
 */
void quicksort_keyval_int_int(int *key, int *val, int start, int end)
{
    int pivot;
    int i, j, k;

    if (start < end)
    {
        k = choose_pivot(start, end);
        swap_int(&key[start], &key[k]);
        swap_int(&val[start], &val[k]);
        pivot = key[start];

        i = start + 1;
        j = end;
        while (i <= j)
        {
            while ((i <= end) && (key[i] <= pivot))
                i++;
            while ((j >= start) && (key[j] > pivot))
                j--;
            if (i < j)
            {
                swap_int(&key[i], &key[j]);
                swap_int(&val[i], &val[j]);
            }
        }

        // swap two elements
        swap_int(&key[start], &key[j]);
        swap_int(&val[start], &val[j]);

        // recursively sort the lesser key
        quicksort_keyval_int_int(key, val, start, j - 1);
        quicksort_keyval_int_int(key, val, j + 1, end);
    }
}


void quicksort_key_int(int *key, int start, int end)
{
    int pivot;
    int i, j, k;

    if (start < end)
    {
        k = choose_pivot(start, end);
        swap_int(&key[start], &key[k]);
        // swap_int(&val[start], &val[k]);
        pivot = key[start];

        i = start + 1;
        j = end;
        while (i <= j)
        {
            while ((i <= end) && (key[i] <= pivot))
                i++;
            while ((j >= start) && (key[j] > pivot))
                j--;
            if (i < j)
            {
                swap_int(&key[i], &key[j]);
                // swap_int(&val[i], &val[j]);
            }
        }

        // swap two elements
        swap_int(&key[start], &key[j]);
        // swap_int(&val[start], &val[j]);
 
        // recursively sort the lesser key
        quicksort_key_int(key,  start, j-1);
        quicksort_key_int(key,  j+1, end);
    }
}


void prefix_sum(int *input, int length)
{
    if (length == 0 || length == 1)
        return;

    // int old_val, new_val;

    // old_val = input[0];
    // input[0] = 0;
    // for (int i = 1; i < length; i++)
    // {
    //     new_val = input[i];
    //     input[i] = old_val + input[i - 1];
    //     old_val = new_val;
    // }

    for (int i = 1; i < length; i++)
    {
        input[i] = input[i - 1] + input[i];
    }
}

// uint64_t getCurrentTimeMilli()
// {
//     struct timeval time;

//     gettimeofday(&time, NULL);
//     return (uint64_t)(time.tv_sec * INT64_C(1000) + time.tv_usec / 1000);
// }

// int fileIsExist(const char *filename)
// {
//     return access(filename, F_OK);
// }

uint64_t doubleToRawBits(double d)
{
    union
    {
        uint64_t i;
        double f;
    } word;
    word.f = d;
    return word.i;
}

uint64_t max_uint64(uint64_t *array, uint64_t length)
{
    if (length < 0)
    {
        printf("max_uint64 legth must l>= 0\n");
        exit(1);
    }
    uint64_t maximum = array[0];
    for (uint64_t i = 0; i < length; i++)
    {
        maximum = maximum < array[i] ? array[i] : maximum;
    }
    return maximum;
}

double max_double(double *array, uint64_t length)
{
    if (length < 0)
    {
        printf("max_double legth must l>= 0\n");
        exit(1);
    }
    double maximum = array[0];
    for (uint64_t i = 0; i < length; i++)
    {
        maximum = maximum < array[i] ? array[i] : maximum;
    }
    return maximum;
}

uint64_t sum_uint64(uint64_t *array, int len)
{
    uint64_t sum = 0;
    for(int i = 0; i< len;i++)
    {
        sum += array[i];
    }
    return sum;
}

#endif