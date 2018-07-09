######  Fortran Complier  ##################################
#
#
#

objects = \
modules.o stopgm.o init.o cv.o wham.o bluemoon.o main.o 

wham.out : $(objects)
	ifort ${FLAGS} -o wham.out $(objects)

$(objects): %.o : %.F90
	ifort ${FLAGS} -c $<

#
#
############################################################
