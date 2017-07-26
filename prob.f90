subroutine prob(xi_l,dxi,n,xi,l,dt,ncut,ks,nsamp,nsteps,str_xi,k,hist,v,fac1)
implicit none

integer,intent(in)::     n,nsamp,ncut,nsteps
integer,intent(inout)::  k
real*8,intent(in)::      xi_l,l,dt,ks,dxi,xi,fac1
real*8,intent(inout),allocatable:: hist(:),v(:)
character(len=12),intent(in)::     str_xi  ! for file reading purposes

integer::                i,j,m,nc,nat,trash,bin
real*8::                 x, dx
real*8::                 tmp,tmp1,tmp2,tmp3,dist,dist1,dist2
character(len=20)::      input, output, test  ! for file reading purposesa
character(len=*), parameter:: prefix = "traj_", dat_sub = ".dat", test_sub = &
                            & ".test", out_sub = ".output"
!real*8, allocatable::    pos(:,:)

input = trim(prefix // adjustl(str_xi)) // dat_sub
output = trim(prefix // adjustl(str_xi)) // out_sub
test = trim(prefix // adjustl(str_xi)) // test_sub
write(*,*) input, output
open(10,file=input,status='old')
open(11,file=output,status='unknown')
open(12,file=test,status='unknown')

! the output normalized probability distribution
!$OMP PARALLEL DO
do i = 1, n               ! initial histogram values to zero
   hist(i) = 0.d0
end do
!$OMP END PARALLEL DO

!nc = 0 ! initialize total number of snapshots to zero
k = 0 ! initialize total number of sampled snapshots to zero
x = (0 - xi_l)/fac1
!x = (xi - xi_l)/fac1
write(*, *) x, xi, xi_l
dx = dxi/fac1

!allocate(pos(3,3))
do m=1,ncut
   read(10,*) 
!   read(10,*) 
!   read(10,*) 
!   read(10,*) 
!   read(10,*) 
!   read(10,*) 
!   read(10,*) 
!   read(10,*) 
!   read(10,*) 
end do

do m=ncut+1,nsteps!/nsamp
     if(mod(m,nsamp).ne.0) cycle
     read(10,*) i,dist 
!     read(10,*) trash,i,tmp,dist 
!     read(10,*)
!     read(10,*)
!     read(10,*) trash,pos(1,1),pos(2,1),pos(3,1)! h
!     read(10,*)
!     read(10,*)
!     read(10,*)
!     read(10,*) trash,pos(1,2),pos(2,2),pos(3,2) ! acceptor o
!     read(10,*) trash,pos(1,3),pos(2,3),pos(3,3) ! donor o

!      useless?

!     dist = 0.d0
!     dist1 = 0.d0
!     dist2 = 0.d0
!     do j = 1,3
!        tmp = (pos(j,1) - pos(j,2)) ! distance between 8 and 4
!        tmp1 = (pos(j,1) - pos(j,3)) ! distance between 9 and 4
!        tmp  = tmp - l*dnint(tmp/l) 
!        tmp1  = tmp1 - l*dnint(tmp1/l) 
!        dist1 = dist1 + tmp*tmp
!        dist2 = dist2 + tmp1*tmp1
!     end do
!     dist = (dsqrt(dist1) - dsqrt(dist2))*fac1
     dist = dist + x
!     dist = dist - xi_l  ! shift distance relative lower bound
     dist = dist/dx
     if(dist.ge.n .or. dist.lt.0)  cycle ! discard points greater than xi_u and smaller than xi_l
     bin = 1 + dint(dist) ! locate bin
     k = k + 1
     hist(bin)=hist(bin) + 1.d0 ! place in appropriate bin
     write(12,*) hist(bin),bin
end do
!30   continue
!     write(*,*)'error in reading input file'
!40   continue
     close(10)
     close(12)

!     compute normalized distribution and biasing potential
write(11,'(a)') '# coordinate     potential     probability'
do i = 1, n
   tmp3    = xi_l + dxi*(dble(i)-1d0) ! reaction coordinate in angstrom
   hist(i) = hist(i)/k ! normalized probality at tmp3
   tmp2 = (tmp3-xi)/fac1 
   v(i) = ks/2*tmp2**2  ! biasing window potential
!   v(i)    = v(i)*fac1**2  ! convert potential to a.u.
   write(11,'(4f12.7)') tmp3,v(i),hist(i)
      end do
close(11)
return
end subroutine
