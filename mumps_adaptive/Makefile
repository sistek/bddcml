# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMPI) $(INCMUMPS) -I../multilevel -I../bddcmumps
# libraries needed
LIBS = $(MUMPSLIBPAR) $(SCALAPACKLIB) $(BLACSLIB) $(BLOPEXLINK) \
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
   ../bddcmumps/module_mumps.o \
   ../bddcmumps/module_sm.o \
   ../multilevel/module_utils.o \
   ../multilevel/module_dd.o

bddcmumps.o: \
   module_bddc.o \
   ../bddcmumps/module_sm.o \
   ../multilevel/module_dd.o \
   ../multilevel/module_adaptivity.o

# building binaries
bddcmumps: \
   bddcmumps.o \
   module_bddc.o \
   ../multilevel/module_dd.o \
   ../multilevel/module_adaptivity.o \
   ../lobpcg_driver/lobpcg_driver.o \
   $(OBJLOBPCGDRIVER) \
   ../multilevel/lobpcg_mvecmult_f.o \
   ../bddcmumps/module_mumps.o \
   ../bddcmumps/module_sm.o \
   ../multilevel/module_utils.o 
	$(FC) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o *.mod bddcmumps $(TESTS)
