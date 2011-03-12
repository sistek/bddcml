# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
INC = $(INCMUMPS) $(INCPARMETIS) $(INCBLOPEX) $(INCMPI)
# libraries needed
LIBS = $(LIBMUMPS) $(LIBPARMETIS) $(LIBSCALAPACK) $(LIBBLACS) \
       $(LIBLAPACK) $(LIBBLAS) $(LIBPORD) $(LIBOTHERS)
# $(LIBMPI)$(LIBPTHREAD)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity

# rules

TESTS = \
   test_module_utils \
   test_module_graph \
   test_module_pp \
   test_module_adaptivity \
   test_module_dd \
   test_module_dd2 \
   test_module_levels \
   test_module_mumps \
   test_module_mumps2 \
   test_module_sm \
   test_module_bddc \
   test_module_bddc2

OBJ = \
   lobpcg_mvecmult_f.o \
   module_krylov.o \
   module_krylov_types_def.o \
   module_adaptivity.o \
   module_dd.o \
   module_levels.o \
   module_mumps.o \
   module_graph.o \
   module_pp.o \
   pget_sub_neighbours.o\
   module_parsing.o \
   module_sm.o \
   module_tecplot.o \
   module_utils.o \
   lobpcg_driver.o \
   bddcml_fortran_interface.o

BDDCLIBDIR=../lib
BDDCLIB=$(BDDCLIBDIR)/libbddcml.a

all: \
   modules tests bddclib

modules: $(OBJ)

tests: $(TESTS)

bddclib: $(BDDCLIB)


# general rule s for compiling objects
%.o : %.f90
	$(MPFC) $(MPFFLAGS) $(INC) -c $<

%.o : %.f
	$(MPFC) $(MPFFLAGS) $(INC) -c $<

%.o : %.c
	$(MPCC) $(MPCFLAGS) $(CDEFS) $(INC) -c $<

# module dependencies
module_bddc.o: \
   module_dd.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o 

module_adaptivity.o: \
   module_dd.o \
   module_utils.o 

module_levels.o: \
   module_adaptivity.o \
   module_dd.o \
   module_pp.o \
   module_graph.o \
   module_parsing.o \
   module_krylov_types_def.o \
   module_utils.o 

module_krylov.o: \
   module_krylov_types_def.o \
   module_levels.o \
   module_utils.o

module_mumps.o: \
   module_utils.o

module_dd.o: \
   module_pp.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o 

module_sm.o: \
   module_utils.o 

module_graph.o: \
   module_parsing.o \
   module_utils.o 

module_pp.o: \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_utils.o 

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

test_module_adaptivity.o: \
   module_pp.o \
   module_adaptivity.o \
   module_dd.o

test_module_levels.o: \
   module_levels.o \
   module_adaptivity.o \
   module_pp.o \
   module_graph.o \
   module_parsing.o \
   module_dd.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o 

test_module_dd.o: \
   module_dd.o

test_module_dd2.o: \
   module_dd.o

test_module_graph.o: \
   module_graph.o

test_module_pp.o: \
   module_pp.o

test_module_utils.o: \
   module_utils.o

lobpcg_mvecmult_f.o: module_levels.o

# building binaries
test_module_mumps: \
   test_module_mumps.o \
   module_mumps.o \
   module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_mumps2: \
   test_module_mumps2.o \
   module_mumps.o \
   module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_sm: \
   test_module_sm.o \
   module_sm.o \
   module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_bddc: \
   test_module_bddc.o \
   module_bddc.o \
   module_dd.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_bddc2: \
   test_module_bddc2.o \
   module_bddc.o \
   module_dd.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_adaptivity: \
   test_module_adaptivity.o \
   module_levels.o \
   module_adaptivity.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_dd.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o \
   module_krylov_types_def.o \
   lobpcg_driver.o \
   lobpcg_mvecmult_f.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS) $(BLOPEXLINK)

test_module_dd: \
   test_module_dd.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_dd.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o 
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_dd2: \
   test_module_dd2.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_dd.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o 
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_levels: \
   test_module_levels.o \
   module_levels.o \
   module_adaptivity.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_dd.o \
   module_mumps.o \
   module_sm.o \
   module_utils.o \
   module_krylov_types_def.o \
   lobpcg_driver.o \
   lobpcg_mvecmult_f.o
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS) $(BLOPEXLINK)

test_module_graph: \
   test_module_graph.o \
   module_graph.o \
   module_parsing.o \
   module_utils.o 
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_pp: \
   test_module_pp.o \
   module_pp.o \
   pget_sub_neighbours.o \
   module_graph.o \
   module_parsing.o \
   module_utils.o 
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

test_module_utils: \
   test_module_utils.o \
   module_utils.o 
	$(MPFC) $(MPFFLAGS) -o $@ $^ $(LIBS)

# make modules into BDDCML archive library
$(BDDCLIB): $(OBJ)
	mkdir -p $(BDDCLIBDIR)
	$(AR) $@ $^
	$(RANLIB) $@

clean:
	-rm -f *.mod *.o $(OBJ) $(TESTS) $(BDDCLIB)
