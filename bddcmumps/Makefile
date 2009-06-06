# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMPI) $(INCMUMPS) -I../common_blocks/
# libraries needed
LIBS = $(LIBMUMPS) $(LIBSCALAPACK) $(LIBBLACS) \
       $(LIBPORD) $(LIBMETIS) $(LIBLAPACK) $(LIBBLAS) $(LIBMPI) \
       $(LIBPTHREAD)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules
all: \
   bddcmumps

# general rule s for compiling objects
%.o : %.f90
	$(MPFC) $(MPFFLAGS) $(INC) -c $<

# module dependencies
module_bddc.o: \
   ../common_blocks/module_mumps.o \
   ../common_blocks/module_sm.o

bddcmumps.o: \
   module_bddc.o \
   ../common_blocks/module_mumps.o \
   ../common_blocks/module_sm.o \
   ../common_blocks/module_utils.o

# building binaries
bddcmumps: \
   bddcmumps.o \
   module_bddc.o \
   ../common_blocks/module_mumps.o \
   ../common_blocks/module_sm.o \
   ../common_blocks/module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o *.mod bddcmumps
