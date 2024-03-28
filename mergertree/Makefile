include ../common.mk

CFLAGS := $(filter-out -Wstrict-aliasing,$(CFLAGS))

SOURCES    = defs.c read_param.c  impulse.c maketree.c loadfillprogenitors.c ../io/loadsnapshot.c ../io/loadgroups.c  ../io/loadparents.c  \
			genplotdata.c ../utils/progressbar.c ../utils/utils.c ../utils/set_cosmology.c main.c

OBJECTS    = $(SOURCES:.c=.o)

EXECUTABLE = MakeMergerTree

INCL       = defs.h read_param.h impulse.h maketree.h loadfillprogenitors.h ../io/io.h ../io/loadsnapshot.h ../io/loadgroups.h ../io/loadparents.h \
						genplotdata.h ../utils/progressbar.h ../utils/utils.h ../utils/set_cosmology.h proto.h 

all: $(SOURCES) $(EXECUTABLE) $(INCL)

$(EXECUTABLE): $(OBJECTS) $(INCL)
	$(CC) $(OBJECTS) $(CLINK) $(GSL_LINK) -o $@

%.o: %.c $(INCL)
	$(CC) $(GSL_CFLAGS) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan celna

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clena: clean
celan: clean
celna: clean


