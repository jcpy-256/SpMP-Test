#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "utils.h"
#include "common_macros.h"

#ifndef MAT_VAL_TYPE
#define MAT_VAL_TYPE double
#endif

#ifndef BENCH_REPEAT
#define BENCH_REPEAT 100
#endif

#ifndef MAT_PTR_TYPE
#define MAT_PTR_TYPE int
#endif


// 纳秒
#ifndef NANOSECOND // 1.8G
#define NANOSECOND (1e9)
#endif

// 微妙
#ifndef MICROSECOND 
#define MICROSECOND (1e6)
#endif

// 毫秒
#ifndef MILLISECOND 
#define MILLISECOND (1e3)
#endif

