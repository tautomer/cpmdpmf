program wham 
!use omp_lib
implicit none

interface
subroutine prob(xi_min,dxi,n,xi,l,dt,ncut,ks,nsamp,nsteps,str_xi,k,hist,v,fac1)
integer,intent(in)::     n,nsamp,nsteps,ncut
integer,intent(inout)::  k
real*8,intent(in)::      l,dt,ks,xi,xi_min,dxi,fac1
real*8,intent(inout),allocatable:: hist(:) ! hist(nmax2) histogram of a given simulation 
real*8,intent(inout),allocatable:: v(:) ! v(nmax2) 
character(len=12),intent(in)::     str_xi ! for file reading purposes
end subroutine prob
end interface

integer::            nw,i,j,k,xi_i
integer::            n,nsamp,ncut,nsteps
real*8::             xi_0, d, xi_max, dxi, kb, temp, beta
real*8::             l,dt,ks, xi_min, xi, fac1
!real*8,allocatable:: ni(:)       ! ni(namx1) smapled snapshots in a given window
real*8,allocatable:: w(:,:)      ! biasing potential
real*8,allocatable:: p_biased(:,:) ! p_biased(namx1,namx2) biased distribution 
real*8,allocatable:: p_unbiased(:) ! p_unbiased(nmax2) 
real*8,allocatable:: fi(:) ! fi(namx1) = exp(-fi*beta) 
real*8,allocatable:: ni(:)
real*8,allocatable:: hist(:) ! hist(nmax2) histogram of a given simulation 
real*8,allocatable:: v(:) ! v(nmax2) 
real,allocatable:: dx(:) 
real*8::             eps, tol, fi_old, pmin, tmp1,tmp2 
real*8::             fac, numerator, denominator, fi_new 
character(len=12)::  str_xi ! convert xi_i to string
character(len=30)::  date

!*************************constants**************************
kb=0.316679d-5       ! boltzmann's constant in a.u./kelvin
fac=627.509d0         ! conversion factor from a.u. to kcal/mol
fac1 = 0.52917721092  ! conversion factor from a.u. to angstrom
!************************************************************

call fdate(date)
open(10,file='free_ener.dat',status='unknown')
write(10, *) '# program started on ', date
!write ( 10, '(a,i8)' ) &
!    '# the number of processors available = ', omp_get_num_procs ()
!write ( 10, '(a,i8)' ) &
!    '# the number of threads available    = ', omp_get_max_threads ()
close(10)

open(unit=20,file="config", action="read")
read(20,*)
read(20,*) xi_0
read(20,*)
read(20,*) xi_min
read(20,*)
read(20,*) xi_max
read(20,*)
read(20,*) dxi
read(20,*)
read(20,*) l
read(20,*)
read(20,*) nsteps
read(20,*)
read(20,*) ncut
read(20,*)
read(20,*) ks
read(20,*)
read(20,*) nsamp
read(20,*)
read(20,*) temp
read(20,*)
read(20,*) nw
read(20,*)
read(20,*) tol
close(20)

n = 1+dint((xi_max-xi_min)/dxi) ! total number of bins
!dt = dt*2.418884d-5 ! convert time from a.u. to ps
beta = 1.d0/kb/temp
!d = (xi_max-xi_min)/(nw)

allocate(dx(nw),ni(nw),p_biased(nw,n),w(nw,n),v(n),hist(n),fi(nw),p_unbiased(n))

dx = 4d-2
!do i = 1, 25
!   dx(i) = 6d-2
!end do
!
!do i = 26, nw
!   dx(i) = 2d-2
!end do

xi = -xi_0
do i = 1, nw ! loop over windows
   xi_i=nint((0.72-xi)*100)
   write(*,*) xi_i
  
   write(str_xi, *) xi_i
   call prob(xi_min,dxi,n,xi,l,dt,ncut,ks,nsamp,nsteps,str_xi,k,hist,v,fac1)
  
   ni(i)=dble(k)
  
!$omp parallel do private(j) shared(p_biased, w)
   do j = 1, n ! loop over historgram points in window i
        p_biased(i,j)=hist(j)  ! biased distribution of window i at coordinate xi_j
        w(i,j)=dexp(-beta*v(j)) ! restraining potential of window i at coordinate xi_j
     end do
!$omp end parallel do
  
   fi(i) = 1.d0 ! initial guess of unity for fi = dexp(-fi*beta)
   xi = xi - dx(i) ! compute central coordinate of window i
end do
eps = tol
do while(eps.ge.tol)
   !$omp parallel do private(j) shared(ni, w, fi, p_biased, p_unbiased) reduction(+:denominator, numerator)
   do j = 1, n
      numerator   = 0.d0
      denominator = 0.d0
      !computing the denominator and numerator
      do i =1,nw
         denominator = denominator + ni(i)*w(i,j)/fi(i)
         numerator = numerator + ni(i)*p_biased(i,j)
      end do
      tmp1=numerator/denominator
      !convert zero probability to very small positive number
      if(tmp1.eq.0.d0) tmp1=1.d-15
      p_unbiased(j)=tmp1
   end do
   !$omp end parallel do
   
   !compute new fi based and the old use eps=sum{(1-fi_new/fi_old)**2}
   eps = 0.d0
   open(13,file='eps.dat',status="unknown")
   do i = 1, nw
      fi_old = fi(i)
      fi_new = 0.d0
      do j = 1, n
         fi_new = fi_new + w(i,j)*p_unbiased(j)
      end do
      fi(i) = fi_new
      eps = eps + (1.d0 - fi_new/fi_old)**2
   end do
   write(13,*) eps
end do
write(*,*) 'converged'
close(13)

!find free energy shift
pmin = -1.d+8
do j = 1, n
   pmin=max(pmin,(1.d0/beta)*dlog(p_unbiased(j)))
end do

!print pmf
open(10,file='free_ener.dat', access='append')
do j = n, 1, -1
   tmp1= xi_min + dxi*(dble(j)-1d0) ! reaction coordinate in angstrom
   write(*,*) tmp1
   !if(tmp1.lt.0) exit
   tmp2=(-(1.d0/beta)*dlog(p_unbiased(j))+pmin)*fac ! pmf in kcal/mol
   write(10, '(2f12.7)') tmp1,tmp2
end do
!do j = 1, n
!   tmp1= xi_min + dxi*(dble(j)-1d0) ! reaction coordinate in angstrom
!   tmp1 = - tmp1
!   if(tmp1.ge.0) cycle
!   tmp2=(-(1.d0/beta)*dlog(p_unbiased(j))+pmin)*fac ! pmf in kcal/mol
!   write(10, '(2f12.7)') tmp1,tmp2
!end do

call fdate(date)
write(10, *) '# program ended on ', date
close(10)
write(*,*) 'pmf printed to free_ener.dat'
end program
