CFLAGS=     # options for C compiler
FFLAGS=
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../hypercube

OPTM =  -g


VPATH = ../01_data_structures:../02_parsers:../03_output:\
	../04_almt_pp:../05_scoring:../06_main: \
        ../07_utils

PLATFORM = $(shell uname)



LIBS    = -lm -lc -llapack
INCLUDE = -I ../00_include


DEFINES =
SOURCES.F =
SOURCES.f =
SOURCES.c = hypercube.c  hc_aa_code.c  hc_afa_out.c  hc_almt.c \
	 hc_comp_score_out.c  hc_chi2.c 	hc_cons_scoring.c \
	 hc_patch_almt.c  hc_rate2prob.c hc_read_alignment.c   hc_read_dssp.c   hc_read_extern.c  \
	 hc_read_groups.c  hc_read_opt_file.c   hc_read_pdb.c \
	 hc_read_tillier.c hc_score_out.c \
	 hc_spec_detection.c  hc_struct_almt_mapping.c hc_tree.c hc_utils.c


AUX =  makefile  hypercube  hc_alignment.h	hc_geometry.h  hc_options.h  hc_pdb.h  hc_tree.h  hc_utils.h  hypercube.h

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM) \
	$(LIBS)    $(INCLUDE)


.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX)
