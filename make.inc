#*********** Compilers
MPFC = gfortran 
MPFFLAGS = -g -C -Wall -fbounds-check
#MPFFLAGS = -O3 
#MPFFLAGS =-g -warn all -check all  
FC = gfortran 
FFLAGS = -g -C -Wall -fbounds-check

MPCC = gcc
MPCFLAGS = -g -Wall
CC = gcc
CFLAGS = -g -Wall

#**********Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#*********** Includes
INCMUMPS = -I/home/sistek/programy/MUMPS_4.8.4/src
INCMPI   = -I/usr/lib/mpich/include

#*********** Libraries
LIBMUMPS     = -L/home/sistek/lib -ldmumps -lmumps_common
LIBSCALAPACK = -L/home/sistek/lib -lscalapack
LIBBLACS     = -L/home/sistek/lib -lblacs -lblacsCinit -lblacsF77init -lblacs 
LIBBLAS      = -L/home/sistek/lib -lblas
LIBLAPACK    = -L/home/sistek/lib -llapack
LIBMETIS     = -L/home/sistek/lib -lmetis
LIBPORD      = -L/home/sistek/lib -lpord
LIBMPI       = -L/usr/lib/mpich/lib -lmpich
LIBMPIFAKE   = -L/home/sistek/lib -lmpiseq
LIBPTHREAD   = -lpthread
LIBOTHERS    = -lm 

# BLOPEX paths
BLOPEXDIR = /home/sistek/programy/lobpcg/blopex
INCBLOPEX = \
   -I$(BLOPEXDIR)/blopex_serial/multivector \
   -I$(BLOPEXDIR)/blopex_serial/pcg_multi \
   -I$(BLOPEXDIR)/blopex_serial/matmultivec \
   -I$(BLOPEXDIR)/blopex_abstract/include
OBJBLOPEX = \
   $(BLOPEXDIR)/blopex_serial/multivector/multi_vector.o \
   $(BLOPEXDIR)/blopex_serial/pcg_multi/pcg_multi.o \
   $(BLOPEXDIR)/blopex_serial/matmultivec/matmultivec.o
LIBBLOPEX = -L$(BLOPEXDIR)/blopex_abstract/lib -lBLOPEX
BLOPEXLINK = $(OBJBLOPEX) $(LIBBLOPEX)

