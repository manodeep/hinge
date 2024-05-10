#include "io.h"

#include "io_hinge_ascii.h"
#include "io_subfind_binary.h"

#ifdef FOF_ONLY
#error "The code does *NOT* work with FOF_ONLY enabled"
#endif

const char GROUP_FORMAT_NAMES[][MAXTAGLEN] = {"subfind_binary", "hinge_ascii"};
const enum valid_group_formats GROUP_FORMAT_ENUMS[] = {subfind_binary, hinge_ascii};
const int nvalid_group_format_names  = sizeof(GROUP_FORMAT_NAMES)/(MAXTAGLEN*sizeof(char));
