######  Fortran Complier  ##################################
#
#
#
FLAGS=-g -O0 -Wall -fcheck=all -fbacktrace

objects = \
modules.o main.o hist.o init.o

wham.out : $(objects)
	gfortran ${FLAGS} -o wham.out $(objects)

$(objects): %.o : %.f90
	gfortran ${FLAGS} -c $<

#
#
############################################################
