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
            perror("System error-msg");                                                                                \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define TEST_FWRITE(ptr, size, nmemb, stream)                                                                          \
    do                                                                                                                 \
    {                                                                                                                  \
        XASSERT(ptr != NULL, "Error: pointer argument for fwrite must not be NULL (requested to write %lld elements of size %zu)\n", (long long) nmemb, size);        \
        XASSERT(stream != NULL, "Error: stream argument for fwrite must not be NULL (requested to write %lld elements of size %zu)\n", (long long) nmemb, size);      \
        XASSERT(size > 0, "Error: size = %zu bytes for fwrite must be greater than 0 (requested to write %lld elements)\n", size, (long long) nmemb);        \
        XASSERT(nmemb > 0, "Error: nmemb argument for fwrite must be greater than 0 (requested to write %lld elements of size %zu)\n", (long long) nmemb, size);      \
        size_t nwritten = fwrite(ptr, size, nmemb, stream);                                                            \
        if (nwritten != (size_t) nmemb)                                                                                \
        {                                                                                                              \
            fprintf(stderr, "Error: fwrite failed to write %lld elements of size %zu to stream\n", (long long) nmemb, size);        \
            perror("System error-msg");                                                                                \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#endif