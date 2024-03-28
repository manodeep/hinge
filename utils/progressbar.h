#ifndef _PROGRESSBAR_H_
#define _PROGRESSBAR_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"

#ifndef MAXLEN
#define MAXLEN 1000
#endif

void init_my_progressbar(const int64_t N, int *interrupted);
void my_progressbar(const int64_t curr_index, int *interrupted);
void finish_myprogressbar(int *interrupted);

#endif
