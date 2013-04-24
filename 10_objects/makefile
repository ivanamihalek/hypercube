CFLAGS=     # options for C compiler
FFLAGS=
FCOMPL= f90   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = ../hyper_c

OPTM =  -O


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
	 hc_patch_almt.c  hc_rate2prob.c hc_read_cw.c   hc_read_dssp.c   hc_read_extern.c  \
	 hc_read_groups.c  hc_read_opt_file.c   hc_read_pdb.c \
	 hc_read_tillier.c hc_score_out.c \
	 hc_spec_detection.c  hc_struct_almt_mapping.c hc_tree.c hc_utils.c




AUX =  makefile  cube  cube_alignment.h  cube_options.h cube_utils.h

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


