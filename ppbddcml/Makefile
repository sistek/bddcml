# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = -I../common_blocks/
# libraries needed
LIBS = $(LIBPARMETIS) $(LIBLAPACK) $(LIBBLAS)
#$(LIBMPI) $(LIBPTHREAD)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules
all: \
   ppbddcml

# general rule s for compiling objects
%.o : %.f90
	$(MPFC) $(MPFFLAGS) $(INC) -c $<

# module dependencies
ppbddcml.o: \
   ../common_blocks/module_pp.o \
   ../common_blocks/module_graph.o \
   ../common_blocks/module_parsing.o \
   ../common_blocks/module_utils.o

# building binaries
ppbddcml: \
   ppbddcml.o \
   ../common_blocks/module_pp.o \
   ../common_blocks/pget_sub_neighbours.o\
   ../common_blocks/module_graph.o \
   ../common_blocks/module_parsing.o \
   ../common_blocks/module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS) 

clean:
	-rm -f *.o *.mod ppbddcml
