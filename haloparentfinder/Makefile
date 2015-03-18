include ../common.mk


EXEC   = HaloParentFinder

SRC    = hierarchy.c  findallparents.c switchfof.c findprogenitor.c read_param.c main.c ../utils/utils.c ../utils/progressbar.c ../io/loadgroups.c 
OBJS  = $(SRC:.c=.o)

INCL   = hierarchy.h findallparents.h switchfof.h findprogenitor.h read_param.h defs.h ../utils/utils.h ../utils/progressbar.h ../io/loadgroups.h ../io/io.h Makefile

all = $(EXEC) $(SRC) $(INCL)

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(CLINK)   -o  $@

$(OBJS): %.o: %.c $(INCL) 
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE)  -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: clean clena celan celna

clena: clean
celan: clean
celna: clean

