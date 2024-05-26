#pragma once

#ifdef NDEBUG
#define XASSERT(EXP, ...)                                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
    } while (0)
#else
#define XASSERT(EXP, ...)                                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        if (!(EXP))                                                                                                    \
        {                                                                                                              \
            fprintf(stderr, "Error in file: %s\tfunc: %s\tline: %d with expression `" #EXP "'\n", __FILE__,            \
                    __FUNCTION__, __LINE__);                                                                           \
            fprintf(stderr, __VA_ARGS__);                                                                              \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)
#endif