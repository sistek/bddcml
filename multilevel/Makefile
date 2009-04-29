# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
#    parallel
 INC = $(INCMPI) $(INCMUMPS) -I../bddcmumps
#    sequential
#INC = $(INCMPIFAKE) $(INCMUMPS)
# libraries needed
#    parallel
#LIBS = $(MUMPSLIBPAR) $(SCALAPACKLIB) $(BLACSLIB) \
       $(LAPACKLIB) $(BLASLIB) $(ORDERINGSLIB) $(LIBMPI) $(LIBOTHERS)
#    semi-sequential
LIBS = $(LIBMPI) $(MUMPSLIBSEQ) $(BLOPEXLINK) \
       $(LAPACKLIB) $(BLASLIB) $(ORDERINGSLIB) $(LIBOTHERS) 
OBJLOBPCGDRIVER = ../lobpcg_driver/lobpcg_driver.o
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

TESTS = \
   test_module_dd \
   test_module_adaptivity \
   test_module_utils 

all: \
   $(TESTS)

# general rule s for compiling objects
%.o : %.f90
	$(FC) $(INC) -c $<

%.o : %.f
	$(FC) $(INC) -c $<

# module dependencies
module_dd.o: \
   module_mumps_seq.o \
   ../bddcmumps/module_sm.o \
   module_utils.o 

module_adaptivity.o: \
   module_dd.o \
   module_utils.o 

lobpcg_mvecmult_f.o: \
   module_dd.o \
   module_utils.o 

test_module_dd.o: \
   module_dd.o

test_module_adaptivity.o: \
   module_adaptivity.o

test_module_utils.o: \
   module_utils.o

# building binaries

$(OBJLOBPCGDRIVER) :
	(cd ../lobpcg_driver ; $(MAKE) )

test_module_dd: \
   test_module_dd.o \
   module_dd.o \
   module_mumps_seq.o \
   ../bddcmumps/module_sm.o \
   module_utils.o
	$(FC) -o $@ $^ $(LIBS)

test_module_adaptivity: \
   test_module_adaptivity.o \
   module_adaptivity.o \
   module_dd.o \
   module_mumps_seq.o \
   ../bddcmumps/module_sm.o \
   module_utils.o \
   $(OBJLOBPCGDRIVER) \
   lobpcg_mvecmult_f.o
	$(FC) -o $@ $^ $(LIBS)

test_module_utils: \
   test_module_utils.o \
   module_utils.o
	$(FC) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o *.mod $(TESTS)
