include ../common.mk

EXECUTABLE := OrphanFixer

SOURCES := maketree.c fillprogenitors.c missinghalos.c read_param.c ../utils/utils.c ../utils/progressbar.c ../utils/utils_read_params.c ../utils/set_cosmology.c
INCL := $(SOURCES:.c=.h)

IO_SRC := ../io/io.c ../io/loadparents.c ../io/loadgroups.c ../io/io_subfind_binary.c ../io/io_hinge_ascii.c ../io/io_hinge_binary.c ../io/io_hinge_utils.c
INCL += $(IO_SRC:.c=.h) ../utils/hinge.h Makefile ../common.mk

SOURCES += $(IO_SRC) main.c
OBJECTS := $(SOURCES:.c=.o)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(CLINK) $(GSL_LINK) -o  $@

$(OBJECTS): $(INCL)
%.o: %.c $(INCL)
	$(CC) $(GSL_CFLAGS) $(INCLUDE) $(CFLAGS) $(OPT) -c $< -o $@

.PHONY: clean clena celan celna

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clena: clean
celan: clean
celna: clean
