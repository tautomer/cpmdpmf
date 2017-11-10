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
sed -i "6s/.*/FLAGS=$FLAGS/" makefile
make
rm -f *.o
