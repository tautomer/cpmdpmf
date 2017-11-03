######  Fortran Complier  ##################################
#
#
#
#FLAGS=-fopenmp

objects = \
modules.o wham.o prob.o readinput.o

wham.out : $(objects)
	gfortran ${FLAGS} -o wham.out $(objects)

$(objects): %.o : %.f90
	gfortran ${FLAGS} -c $<

#
#
############################################################
