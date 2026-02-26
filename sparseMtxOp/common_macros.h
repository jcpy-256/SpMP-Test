#ifndef COMMON_MACROS_H
#define COMMON_MACROS_H

#define FREE_PTR(ptr)      \
    do                     \
    {                      \
        if ((ptr) != NULL) \
        {                  \
            free(ptr);     \
            (ptr) = NULL;  \
        }                  \
    } while (0)


    
#define CHECK__PTR(ptr)                                        \
    do                                                         \
    {                                                          \
        if ((ptr) == NULL)                                     \
        {                                                      \
            fprintf(stderr, "[ERROR] Null pointer detected!\n" \
                            "  File: %s\n"                     \
                            "  Line: %d\n"                     \
                            "  Function: %s\n",                \
                    __FILE__, __LINE__, __func__);             \
            exit(EXIT_FAILURE); /* 直接退出程序 */             \
        }                                                      \
    } while (0)

#endif