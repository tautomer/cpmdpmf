######  Fortran Complier  ##################################
#
#
#

objects = \
modules.o main.o hist.o init.o

wham.out : $(objects)
	ifort ${FLAGS} -o wham.out $(objects)

$(objects): %.o : %.f90
	ifort ${FLAGS} -c $<

#
#
############################################################
