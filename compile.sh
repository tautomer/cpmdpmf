flag1='-fopenmp'
flag2='-g -check all -fpe0 -warn -traceback -debug extended'
if [[ $1 -eq 1 ]] 
then
	FLAGS=$flag1
elif [[ $1 -eq 2 ]] 
then
	FLAGS=''
elif [[ $1 -eq 3 ]] 
then
	FLAGS="$flag1 $flag2"
else
	FLAGS=$flag2
fi
make FLAGS=$FLAGS
rm -f *.o
