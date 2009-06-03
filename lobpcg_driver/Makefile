# general setting on system (compiler, flags, libraries)
include ../make.inc

# libraries needed

INC=$(BLOPEXINC)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

all: lobpcg_driver.o

# general rule s for compiling objects
%.o : %.f90
	$(FC) -c $<

%.o : %.c
	$(CC) $(INC) $(CDEFS) -c $<

clean:
	rm -f *~ *.o 
