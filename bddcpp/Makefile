# general setting on system (compiler, flags, libraries)
include ../make.inc

# libraries needed
LIBS = $(ORDERINGSLIB)
# includes
INC = -I../bddcmumps
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

all: bddcpp bddcpp2

# general rule s for compiling objects
%.o : %.f90
	$(FC) $(INC) -c $<

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

bddcpp2: \
   bddcpp2.o \
   ../bddcmumps/module_utils.o \
   module_parsing.o \
   module_tecplot.o \
   module_errors.o 
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f *~ *.o *.mod bddcpp bddcpp2
