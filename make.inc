#*********** Compiler
#FC = gfortran -O3 -fbounds-check
FC = gfortran -g -C -Wall -fbounds-check

CC = gcc -g -Wall

#FC=ifort -g -warn all -check all  

#FC=g95 -g -Wall

#*********** Includes
INCMUMPS = -I/home/sistek/programy/MUMPS_4.8.4/src
INCMPI = -I/usr/lib/mpich/include
INCMPIFAKE = -I/home/sistek/programy/MUMPS_4.8.4/libseq

#*********** Libraries
MUMPSLIBPAR = -L/home/sistek/lib -ldmumps -lmumps_common
MUMPSLIBSEQ = -L/home/sistek/lib -ldmumpsseq 
SCALAPACKLIB = -L/home/sistek/lib -lscalapack
BLACSLIB = -L/home/sistek/lib -lblacs -lblacsCinit -lblacsF77init -lblacs 
BLASLIB = -L/home/sistek/lib -lblas
LAPACKLIB = -L/home/sistek/lib -llapack
ORDERINGSLIB = -L/home/sistek/lib -lmetis -lpord
LIBMPI = -L/usr/lib/mpich/lib -lmpich
LIBMPIFAKE = -L/home/sistek/lib -lmpiseq
LIBOTHERS = -lpthread

 #**********Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

BLOPEX=/home/sistek/programy/lobpcg/blopex
BLOPEXLIB=-L$(BLOPEX)/blopex_abstract/lib -lBLOPEX
BLOPEXOBJ=$(BLOPEX)/blopex_serial/multivector/multi_vector.o $(BLOPEX)/blopex_serial/pcg_multi/pcg_multi.o $(BLOPEX)/blopex_serial/matmultivec/matmultivec.o
BLOPEXINC=-I$(BLOPEX)/blopex_serial/multivector -I$(BLOPEX)/blopex_serial/pcg_multi -I$(BLOPEX)/blopex_serial/matmultivec -I$(BLOPEX)/blopex_abstract/include
LIBC=-lm 
BLOPEXLINK=$(BLOPEXOBJ) $(BLOPEXLIB)

