# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
 INC = $(INCMPI) $(INCMUMPS) -I../multilevel/
# libraries needed
LIBS = $(MUMPSLIBPAR) $(SCALAPACKLIB) $(BLACSLIB) \
       $(LAPACKLIB) $(BLASLIB) $(ORDERINGSLIB) $(LIBMPI) $(LIBOTHERS)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

TESTS = \
   test_module_mumps \
   test_module_mumps2 \
   test_module_sm \
   test_module_bddc \
   test_module_bddc2

all: \
   bddcmumps \
   $(TESTS)

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

test_module_mumps.o: \
   module_mumps.o

test_module_mumps2.o: \
   module_mumps.o

test_module_sm.o: \
   module_sm.o

test_module_bddc.o: \
   module_bddc.o

test_module_bddc2.o: \
   module_bddc.o

# building binaries
bddcmumps: \
   bddcmumps.o \
   module_bddc.o \
   module_mumps.o \
   module_sm.o \
   condest.o
	$(FC) -o $@ $^ $(LIBS)

test_module_mumps: \
   test_module_mumps.o \
   module_mumps.o
	$(FC) -o $@ $^ $(LIBS)

test_module_mumps2: \
   test_module_mumps2.o \
   module_mumps.o
	$(FC) -o $@ $^ $(LIBS)

test_module_sm: \
   test_module_sm.o \
   module_sm.o
	$(FC) -o $@ $^ $(LIBS)

test_module_bddc: \
   test_module_bddc.o \
   module_bddc.o \
   module_mumps.o \
   module_sm.o
	$(FC) -o $@ $^ $(LIBS)

test_module_bddc2: \
   test_module_bddc2.o \
   module_bddc.o \
   module_mumps.o \
   module_sm.o
	$(FC) -o $@ $^ $(LIBS)

clean:
	-rm -f *.o *.mod bddcmumps $(TESTS)
