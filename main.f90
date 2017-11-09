program wham
    use global
    implicit none
    !use omp_lib

    integer i, j, k
    ! w(:,:)          biasing potential
    ! p_biased(:,:)   biased distribution
    ! p_unbiased(:)
    ! fi(:)           fi = exp(-fi*beta)
    ! hist(:)
    ! v(:)
    real*8, allocatable :: hist(:), v(:), w(:,:), p_biased(:,:)
    character(len=30) :: date, cwd
    character(len=100) :: rootdir

    call fdate(date)
    call read_conf()
    call read_meta()
    call init_param(date)
    allocate(p_biased(nw , n), w(nw, n), v(n), hist(n))
    call getcwd(rootdir)
    write(*, *) rootdir

    do i = 1, nw ! loop over windows
        cwd = trim("./" // adjustl(dir(i)))
        call chdir(cwd)
        call prob(i, hist, v)
        call chdir(rootdir)
        ! biased distribution of window i at coordinate xi_j
        p_biased(i, :) = hist
        ! restraining potential of window i at coordinate xi_j
        w(i, :) = dexp(-beta * v)
    end do

    call unbias(w, p_biased)

end program

subroutine unbias(w, p_biased)
    use global
    implicit none

    integer i, j, k
    real*8 eps, fi_old, pmin, tmp1, tmp2
    real*8 numerator, denominator, fi_new, invbeta
    real*8 p_unbiased(n), fi(nw), tmp(n)
    character(len=30) ::  date
    real*8, intent(in) :: w(nw, n), p_biased(nw, n)

    fi = 1.d0 ! initial guess of unity for fi = dexp(-fi*beta)
    eps = tol
    k = 0
    do while(eps.ge.tol)
        !$omp parallel do private(j) shared(ni, w, fi, p_biased, p_unbiased) reduction(+:denominator, numerator)
        do j = 1, n
            numerator = 0.d0
            denominator = 0.d0
            !computing the denominator and numerator
            do i =1, nw
                denominator = denominator + ni(i)*w(i,j)/fi(i)
                numerator = numerator + ni(i)*p_biased(i,j)
            end do
            tmp1 = numerator / denominator
            !convert zero probability to very small positive number
            if(tmp1.eq.0) tmp1 = 1.d-15
            p_unbiased(j) = tmp1
        end do
        !$omp end parallel do

        !compute new fi based and the old use eps=sum{(1-fi_new/fi_old)**2}
        eps = 0.d0
        do i = 1, nw
            fi_old = fi(i)
            fi_new = 0.d0
            do j = 1, n
                fi_new = fi_new + w(i, j) * p_unbiased(j)
            end do
            fi(i) = fi_new
            eps = eps + (1.d0 - fi_new / fi_old) ** 2
        end do
        k = k + 1
        if(mod(k, 10000).eq.1) write(*, *) k, eps
    end do
    write(*,*) 'converged after ', k, ' loops'

    !find free energy shift
    invbeta = 1.d0 / beta
    k = maxloc(p_unbiased, dim=1)
    j = n / 2
    pmin = invbeta * dlog(maxval(p_unbiased))
    if(k.gt.j) then
        tmp = [p_unbiased(n:n-j+1:-1), p_unbiased(j+1:n)]
    else
        tmp = [p_unbiased(1:j), p_unbiased(j+1:1:-1)]
    end if

    !print pmf
    open(10,file='free_ener.dat', access='append')
    do j = 1, n
        tmp1 = xbin(j) / fac
        tmp2 = (-invbeta * dlog(p_unbiased(j)) + pmin) * fac1
        write(10, '(2f12.7)') tmp1, tmp2
    end do

    call fdate(date)
    write(10, *) '# program ended on ', date
    close(10)
    write(*, *) 'pmf printed to free_ener.dat'
end subroutine