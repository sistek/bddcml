# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMPI) $(INCMUMPS) -I../multilevel/
# libraries needed
LIBS = $(MUMPSLIBPAR) $(MUMPSLIBSEQ) $(SCALAPACKLIB) $(BLACSLIB) $(BLOPEXLINK) \
       $(LAPACKLIB) $(BLASLIB) $(ORDERINGSLIB) $(LIBMPI) $(LIBOTHERS)
OBJLOBPCGDRIVER = ../lobpcg_driver/lobpcg_driver.o
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

TESTS = 

all: \
   bddcmumps 

# general rule s for compiling objects
%.o : %.f90
	$(FC) $(INC) -c $<

%.o : %.f
	$(FC) $(INC) -c $<

# module dependencies
module_bddc.o: \
   module_mumps.o \
   module_sm.o

bddcmumps.o: \
   module_bddc.o \
   module_sm.o

# building binaries
bddcmumps: \
   bddcmumps.o \
   module_bddc.o \
   ../multilevel/module_dd.o \
   ../multilevel/module_adaptivity.o \
   ../lobpcg_driver/lobpcg_driver.o \
   $(OBJLOBPCGDRIVER) \
   ../multilevel/lobpcg_mvecmult_f.o \
   module_mumps.o \
   ../multilevel/module_mumps_seq.o \
   module_sm.o \
   condest.o \
   ../multilevel/module_utils.o 
	$(FC) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o *.mod bddcmumps $(TESTS)
