include ../common.mk


OPTIONS =  $(OPTIMIZE) $(OPT) $(CCFLAGS)

SOURCES    = defs.c  maketree.c fillprogenitors.c missinghalos.c read_param.c ../utils/utils.c ../utils/progressbar.c ../utils/set_cosmology.c ../io/loadsnapshot.c ../io/loadgroups.c  ../io/loadparents.c main.c 

OBJECTS    = $(SOURCES:.c=.o)

EXECUTABLE = OrphanFixer

INCL       = defs.h maketree.h fillprogenitors.h missinghalos.h read_param.h ../utils/utils.h ../utils/progressbar.h ../utils/set_cosmology.h ../io/loadsnapshot.h ../utils/sglib.h ../io/io.h ../io/loadgroups.h ../io/loadparents.h Makefile

all: $(SOURCES) $(EXECUTABLE) $(INCL)

$(EXECUTABLE): $(OBJECTS) $(INCL)
	$(CC) $(CLINK) $(GSL_LINK) $(OBJECTS) -o $@

$(OBJECTS): %.o: %.c $(INCL)
	$(CC) $(GSL_CFLAGS) $(CFLAGS) $(OPT) $(INCLUDE) -c $< -o $@

.PHONY: clean clena celan celna

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

clena: clean
celan: clean
celna: clean

