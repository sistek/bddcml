# general setting on system (compiler, flags, libraries)
include ../../make.inc

# libraries needed
LIBS = $(ORDERINGSLIB)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

all: bddcpp

# general rule s for compiling objects
%.o : %.f90
	$(FC) -c $<

# module dependencies
bddcpp.o: \
   module_parsing.o \
   module_tecplot.o \
   module_errors.o 

# building binaries
bddcpp: \
   bddcpp.o \
   module_parsing.o \
   module_tecplot.o \
   module_errors.o 
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f *~ *.o *.mod

veryclean: clean
	rm -f bddcpp
