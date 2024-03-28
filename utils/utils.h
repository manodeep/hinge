#pragma once

#include <assert.h>
#include <inttypes.h> //defines PRId64 for printing int64_t
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h> //defines int64_t datatype -> *exactly* 8 bytes int
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#ifndef DM_PART_TYPE
#define DM_PART_TYPE 1
#endif

#include "io.h"

int my_snprintf(char *buffer, int len, const char *format, ...)
    __attribute__((format(printf, 3, 4)));
int64 getnumlines(const char *fname, const char comment);
float periodic(float dx);
float periodic_wrap(float x);
void print_time(const time_t t0, const time_t t1, const char *s);
FILE *my_fopen(const char *fname, const char *mode);
FILE *my_fopen_carefully(const char *fname, void (*header)(FILE *));
void move_existing_file(const char *fname);
void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void *my_malloc(size_t size, int64 N);
void *my_realloc(void *x, size_t size, int64_t N, const char *varname);

void *my_calloc(size_t size, int64 N);
void my_free(void **x);

// gadget snapshot related utils
int get_gadget_nfiles(const char *fname);
int64 get_Numpart(struct io_header *header);
struct io_header get_gadget_header(const char *fname);

// group struct related utils
void group_init(struct group_data *g, int64 N);
double compute_rank(int64 i);
void free_group(struct group_data *g, int64 N);
void free_group_positions(struct group_data *g, int64 N);
struct group_data *allocate_group(int64 N);

// utils related to the particle matching
void init_all_ranks(double *rank, int64 *ncommon, int64 N);
int64 find_max_rank(double *NextAllRanks, int64 NextNsub);
void reset_ncommon(int64 *Ncommon, int64 N);
int64 find_max_ncommon(int64 *Ncommon, int64 N);
int64 get_ncommon(struct group_data *prev, struct group_data *next);

short float_almost_equal(float A, float B, int maxUlps);
