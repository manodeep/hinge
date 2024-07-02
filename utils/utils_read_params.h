#include "hinge.h"
#include "io.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define INT64 4
#define MAXTAGS 20

extern void read_params(const char *fname, struct params_data *params,
                        void special_params(struct params_data *params, const int maxtags, void **addr, int *id,
                                            char (*tag)[MAXLEN], int *nt));
extern void output_params(const char *fname, struct params_data *params,
                          void output_special_params(FILE *, struct params_data *));
extern void set_simulation_params(struct params_data *params);
