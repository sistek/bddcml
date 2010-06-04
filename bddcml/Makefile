# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMUMPS) -I../common_blocks/
# libraries needed
LIBS = $(LIBMUMPS) $(LIBPARMETIS) $(LIBSCALAPACK) $(LIBBLACS) \
       $(LIBPORD) $(LIBLAPACK) $(LIBBLAS) 
#$(LIBMPI) $(LIBPTHREAD)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules
all: \
   bddcml

# general rule s for compiling objects
%.o : %.f90
	$(MPFC) $(MPFFLAGS) $(INC) -c $<

# module dependencies
bddcml.o: \
   ../common_blocks/module_krylov_types_def.o \
   ../common_blocks/module_levels.o \
   ../common_blocks/module_dd.o \
   ../common_blocks/module_mumps.o \
   ../common_blocks/module_sm.o \
   ../common_blocks/module_utils.o

# building binaries
bddcml: \
   bddcml.o \
   ../common_blocks/module_krylov_types_def.o \
   ../common_blocks/module_levels.o \
   ../common_blocks/module_dd.o \
   ../common_blocks/module_mumps.o \
   ../common_blocks/module_sm.o \
   ../common_blocks/module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o *.mod bddcml
