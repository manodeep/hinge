include ../common.mk


EXEC   = HaloParentFinder

SRC := hierarchy.c  findallparents.c switchfof.c findprogenitor.c read_param.c ../utils/utils.c ../utils/progressbar.c ../utils/utils_read_params.c ../utils/set_cosmology.c
INCL := $(SRC:.c=.h)

IO_SRC := ../io/io.c ../io/loadgroups.c ../io/io_subfind_binary.c ../io/io_hinge_ascii.c ../io/io_hinge_binary.c ../io/io_hinge_utils.c

SRC  += $(IO_SRC) main.c
INCL += $(IO_SRC:.c=.h) ../utils/hinge.h Makefile

OBJS   := $(SRC:.c=.o)

all = $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLINK) $(GSL_LINK) -o  $@

$(OBJS): $(INCL)
%.o: %.c
	@echo "Compiling $<"
	$(CC) $(GSL_CFLAGS) $(INCLUDE) $(CFLAGS) $(OPT) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: clean clena celan celna

clena: clean
celan: clean
celna: clean
