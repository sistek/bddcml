# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMPI) $(INCMUMPS) -I../common_blocks
# libraries needed
LIBS = $(LIBMUMPS) $(LIBSCALAPACK) $(LIBBLACS) \
       $(LIBLAPACK) $(LIBBLAS) $(LIBPORD) $(LIBMETIS) \
       $(LIBMPI) $(LIBPTHREAD)
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
bddcmumps.o: \
   ../common_blocks/module_bddc.o \
   ../common_blocks/module_sm.o \
   ../common_blocks/module_dd.o \
   ../common_blocks/module_adaptivity.o

# building binaries
bddcmumps: \
   bddcmumps.o \
   ../common_blocks/module_bddc.o \
   ../common_blocks/module_dd.o \
   ../common_blocks/module_adaptivity.o \
   ../lobpcg_driver/lobpcg_driver.o \
   ../common_blocks/lobpcg_mvecmult_f.o \
   ../common_blocks/module_mumps.o \
   ../common_blocks/module_sm.o \
   ../common_blocks/module_utils.o 
	$(FC) -o $@ $^ $(LIBS) $(BLOPEXLINK)

clean:
	-rm -f *.o *.mod bddcmumps
