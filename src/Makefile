# general setting on system (compiler, flags, libraries)
include ../make.inc

# include files needed
INC = $(BDDCML_INCMUMPS) $(BDDCML_INCPARMETIS) $(BDDCML_INCMETIS) $(BDDCML_INCBLOPEX) $(BDDCML_INCMPI)
##############################################################
# GNU Make automatic variables:
# $@ - target
# $^ - all prerequisities
# $< - first prerequisity


OBJ = \
   lobpcg_mvecmult_f.o \
   module_bddcml_c_wrapper.o \
   module_bddcml.o \
   module_krylov.o \
   module_krylov_types_def.o \
   module_adaptivity.o \
   module_dd.o \
   module_densela.o \
   module_levels.o \
   module_mumps.o \
   module_graph.o \
   module_pp.o \
   pget_sub_neighbours.o \
   pdivide_mesh.o \
   graph_divide.o \
   module_parsing.o \
   module_sm.o \
   module_paraview.o \
   module_utils.o \
   lobpcg_driver.o 

# rules
all: \
   modules bddclib

modules: $(OBJ)

bddclib: $(LIBBDDCML)

# general rule for compiling objects
%.o : %.f90
	$(BDDCML_MPFC) $(BDDCML_MPFFLAGS) $(INC) -c $<

%.o : %.c
	$(BDDCML_MPCC) $(BDDCML_MPCFLAGS) $(BDDCML_CDEFS) $(INC) -c $<

# module dependencies
module_bddcml_c_wrapper.o: \
   module_bddcml.o

module_bddcml.o: \
   module_levels.o \
   module_krylov.o \
   module_utils.o 

module_adaptivity.o: \
   module_dd.o \
   module_graph.o \
   module_utils.o 

module_levels.o: \
   module_adaptivity.o \
   module_dd.o \
   module_pp.o \
   module_graph.o \
   module_parsing.o \
   module_krylov_types_def.o \
   module_paraview.o \
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
   module_paraview.o \
   module_utils.o 

module_sm.o: \
   module_utils.o 

module_graph.o: \
   module_parsing.o \
   graph_divide.o \
   module_utils.o 

module_pp.o: \
   pget_sub_neighbours.o \
   pdivide_mesh.o \
   module_graph.o \
   module_parsing.o \
   module_utils.o 

module_paraview.o: \
   module_utils.o 

lobpcg_mvecmult_f.o: \
   module_levels.o

# make modules and join them into BDDCML archive library
$(LIBBDDCML): $(OBJ)
	mkdir -p $(BDDCML_LIB_DIR)
	$(BDDCML_AR) $@ $^
	$(BDDCML_RANLIB) $@

clean:
	-rm -f *.mod *.o $(OBJ) $(LIBBDDCML)
