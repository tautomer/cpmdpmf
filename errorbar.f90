program errorbar
implicit none

integer:: i, k, nlines
real:: mean, sd, cv, temp
real,allocatable:: prob(:)

! initialization
nlines = 12301
allocate(prob(5))

! open I/O files
open(11,file='FREE_ENER_1.dat')
open(12,file='FREE_ENER_2.dat')
open(13,file='FREE_ENER_3.dat')
open(14,file='FREE_ENER_4.dat')
open(15,file='FREE_ENER_5.dat')
open(20,file='FREE_ENER.dat')

! loop over steps
do i=1,nlines
! need to be zeroed every step
   mean = 0
   sd = 0
   temp = 0
! loop over 5 files
   do k=1,5
   read(10+k,*) cv, prob(k)
! calculate mean value
   mean = mean + prob(k)
   end do
   mean = mean/5
! naive way to get stanard deviation
   do k=1,5
!  to get squared deviation
   temp = temp + (prob(k) - mean)**2
   end do
!  stanard deviation
!  may need improvement
!  commonly error bar isn't like this right?
   sd = sqrt(temp)/5
   write(20,*) cv,prob,mean,sd
end do
end program
