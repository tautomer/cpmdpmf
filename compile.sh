if [[ $1 -eq 1 ]]
then
    FC='ifort'
	flag1='-qopenmp'
	flag2='-g -check all -fpe0 -warn -traceback -debug extended'
else
    FC='gfortran'
	flag1='-fopenmp'
	flag2='-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none
	       -fbacktrace -ffree-line-length-0 -fcheck=all
	   	   -ffpe-trap=zero,overflow,underflow -finit-real=nan'
fi
if [[ $2 -eq 1 ]]
then
    FLAGS=$flag1
elif [[ $2 -eq 2 ]]
then
    FLAGS=''
elif [[ $2 -eq 3 ]]
then
    FLAGS="$flag1 $flag2"
else
    FLAGS=$flag2
fi
make FC="$FC" FLAGS="$FLAGS"
rm -f *.o *genmod*
