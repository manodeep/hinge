include ../common.mk


EXEC   = HaloParentFinder

IO_SRC := ../io/io.c ../io/loadgroups.c ../io/io_subfind_binary.c ../io/io_hinge_ascii.c ../io/io_hinge_binary.c ../io/io_hinge_utils.c
SRC    := hierarchy.c  findallparents.c switchfof.c findprogenitor.c read_param.c main.c ../utils/utils.c ../utils/progressbar.c ../utils/utils_read_params.c
SRC    += $(IO_SRC)
OBJS   := $(SRC:.c=.o)

INCL   := $(IO_SRC:.c=.h)
INCL   += Makefile defs.h

all = $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(CLINK) -o  $@

$(OBJS): $(INCL)
%.o: %.c
	@echo "Compiling $<"
	$(CC) $(INCLUDE) $(CFLAGS) $(OPT) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: clean clena celan celna

clena: clean
celan: clean
celna: clean
