include ../common.mk

EXECUTABLE = MakeMergerTree

CFLAGS := $(filter-out -Wstrict-aliasing,$(CFLAGS))

SOURCES := read_param.c  impulse.c maketree.c loadfillprogenitors.c genplotdata.c ../utils/progressbar.c ../utils/utils.c ../utils/set_cosmology.c ../utils/utils_read_params.c
INCL := $(SOURCES:.c=.h)

IO_SRC := ../io/io.c ../io/loadparents.c ../io/loadgroups.c ../io/io_subfind_binary.c ../io/io_hinge_ascii.c ../io/io_hinge_binary.c ../io/io_hinge_utils.c
INCL += $(IO_SRC:.c=.h) ../utils/hinge.h Makefile

SOURCES += $(IO_SRC) main.c
OBJECTS := $(SOURCES:.c=.o)

OBJECTS    = $(SOURCES:.c=.o)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(CLINK) $(GSL_LINK) -o $@

$(OBJECTS): $(INCL)
%.o: %.c $(INCL)
	$(CC) $(INCLUDE) $(GSL_CFLAGS) $(CFLAGS) $(OPT)  -c $< -o $@

.PHONY: clean clena celan celna

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clena: clean
celan: clean
celna: clean
