program wham
    use global
    implicit none
    !use omp_lib


    integer            i, j, k
    real*8,allocatable:: w(:,:)      ! biasing potential
    real*8,allocatable:: p_biased(:,:) ! biased distribution
    real*8,allocatable:: p_unbiased(:)
    real*8,allocatable:: fi(:) ! fi = exp(-fi*beta)
    real*8,allocatable:: hist(:)
    real*8,allocatable:: v(:)
    real*8::             eps, fi_old, pmin, tmp1,tmp2
    real*8::             numerator, denominator, fi_new
    character(len=30)::  date


    call fdate(date)
    open(10,file='free_ener.dat',status='unknown')
    write(10, *) '# program started on ', date
    !write ( 10, '(a,i8)' ) &
    !    '# the number of processors available = ', omp_get_num_procs ()
    !write ( 10, '(a,i8)' ) &
    !    '# the number of threads available    = ', omp_get_max_threads ()
    close(10)

    call read_conf()
    allocate(dir(nw), xi(nw), nsteps(nw), ni(nw), fi(nw))
    call read_meta()
    allocate(xbin(n), p_biased(nw , n), w(nw, n), v(n), hist(n), p_unbiased(n))
    !interface
    !    subroutine prob(i, hist, v)
    !        real*8,intent(out) :: hist(n) ! histogram  a given simulation
    !        real*8,intent(out) :: v(n)
    !    end subroutine prob
    !end interface
    ni = 0
    do i = 1, n
        xbin(i) = xmin + wbin * (i - 0.5)
    end do

    do i = 1, nw ! loop over windows
        call prob(i, hist, v)

        !$omp parallel do private(j) shared(p_biased, w)
        do j = 1, n ! loop over historgram points in window i
            ! biased distribution of window i at coordinate xi_j
            p_biased(i, j) = hist(j)
            ! restraining potential of window i at coordinate xi_j
            w(i, j) = dexp(-beta * v(j))
        end do
        !$omp end parallel do

    end do

    fi = 1.d0 ! initial guess of unity for fi = dexp(-fi*beta)
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
    do j = 1, n
        tmp1= xbin(j) / fac
        write(*,*) tmp1
        !if(tmp1.lt.0) exit
        tmp2=(-(1.d0/beta)*dlog(p_unbiased(j))+pmin)*fac1 ! pmf in kcal/mol
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
