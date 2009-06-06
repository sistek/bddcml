# general setting on system (compiler, flags, libraries)
include ../make.inc

# libraries needed

INC=$(INCBLOPEX)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

all: lobpcg_driver.o

# general rule s for compiling objects
%.o : %.c
	$(MPCC) $(MPCFLAGS) $(CDEFS) $(INC) -c $<

clean:
	rm -f *~ *.o 
