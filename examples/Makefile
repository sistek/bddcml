# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMUMPS) $(INCMPI) -I../src 
# libraries needed
LIBS = $(LIBBDDCML) $(LIBMUMPS) $(LIBPARMETIS) $(LIBSCALAPACK) $(LIBBLACS) \
       $(LIBPORD) $(LIBLAPACK) $(LIBBLAS) $(BLOPEXLINK) $(LIBOTHERS)
EXAMPLES = bddcml_global bddcml_local
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity
# rules
all: $(EXAMPLES)

# general rule s for compiling objects
%.o : %.f90
	$(MPFC) $(MPFFLAGS) $(INC) -c $<

# building binaries
bddcml_global: \
   bddcml_global.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBBDDCML) $(LIBS) 

bddcml_local: \
   bddcml_local.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBBDDCML) $(LIBS) 
clean:
	-rm -f *.o *.mod $(EXAMPLES)
