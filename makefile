######  Fortran Complier  ##################################
#
#
#

objects = \
modules.o init.o cv.o wham.o bluemoon.o main.o 

wham.out : $(objects)
	${FC} ${FLAGS} -o wham.out $(objects)

$(objects): %.o : %.F90
	${FC} ${FLAGS} -c $<

#
#
############################################################
