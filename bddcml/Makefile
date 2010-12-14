# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMUMPS) -I../common_blocks/
# libraries needed
LIBS = $(LIBMUMPS) $(LIBPARMETIS) $(LIBSCALAPACK) $(LIBBLACS) \
       $(LIBPORD) $(LIBLAPACK) $(LIBBLAS) $(BLOPEXLINK)
LIBBDDCML = -L../common_blocks -lbddcml
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
   ../common_blocks/libbddcml.a

# building binaries
bddcml: \
   bddcml.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBBDDCML) $(LIBS) 

clean:
	-rm -f *.o *.mod bddcml
