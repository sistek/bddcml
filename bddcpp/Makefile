# general setting on system (compiler, flags, libraries)
include ../make.inc

# libraries needed
LIBS = $(LIBMETIS) $(LIBLAPACK) $(LIBBLAS)
# includes
INC = -I../common_blocks
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

all: bddcpp bddcpp2

# general rule s for compiling objects
%.o : %.f90
	$(FC) $(FFLAGS) $(INC) -c $<

# module dependencies
bddcpp.o: \
   ../common_blocks/module_parsing.o \
   ../common_blocks/module_tecplot.o \
   ../common_blocks/module_utils.o

bddcpp2.o: \
   ../common_blocks/module_parsing.o \
   ../common_blocks/module_tecplot.o \
   ../common_blocks/module_utils.o

# building binaries
bddcpp: \
   bddcpp.o \
   ../common_blocks/module_parsing.o \
   ../common_blocks/module_tecplot.o \
   ../common_blocks/module_utils.o 
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

bddcpp2: \
   bddcpp2.o \
   ../common_blocks/module_parsing.o \
   ../common_blocks/module_tecplot.o \
   ../common_blocks/module_utils.o 
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *~ *.o *.mod bddcpp bddcpp2
