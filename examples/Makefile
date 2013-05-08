# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
# libraries needed
LIBS = $(BDDCML_LINK_SEQUENCE) 
EXAMPLES = bddcml_global bddcml_local poisson_on_cube
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity
# rules
all: $(EXAMPLES)

# general rule s for compiling objects
%.o : %.f90
	$(BDDCML_MPFC) $(BDDCML_MPFFLAGS) -I../src -c $<

# building binaries
poisson_on_cube: \
   poisson_on_cube.o \
   $(LIBBDDCML)
	$(BDDCML_MPFC) $(BDDCML_MPFFLAGS) -o $@ $^ $(LIBS) 

bddcml_global: \
   bddcml_global.o \
   $(LIBBDDCML)
	$(BDDCML_MPFC) $(BDDCML_MPFFLAGS) -o $@ $^ $(LIBS) 

bddcml_local: \
   bddcml_local.o \
   $(LIBBDDCML)
	$(BDDCML_MPFC) $(BDDCML_MPFFLAGS) -o $@ $^ $(LIBS) 
clean:
	-rm -f *.o *.mod $(EXAMPLES)
