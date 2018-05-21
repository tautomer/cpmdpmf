program wham
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    implicit none

    integer i
    ! w(:,:)          biasing potential
    ! p_biased(:,:)   biased distribution
    ! p_unbiased(:)
    ! fi(:)           fi = exp(-fi*beta)
    ! hist(:)
    ! v(:)
    real*8, allocatable :: hist(:), v(:), w(:,:), p_biased(:,:)
    character(len=30) :: date

    call fdate(date)
    call read_conf()
    call folderloop()
    call init_param(date)
    call movefile(nw, 1)
    allocate(p_biased(nw , n), w(nw, n), v(n), hist(n))

#if defined(_OPENMP)
    !$omp parallel do &
    !$omp private(i, hist, v) &
    !$omp shared(p_biased, w)
#endif
    do i = 1, nw ! loop over windows
        call prob(i, hist, v)
        ! biased distribution of window i at coordinate xi_j
        p_biased(i, :) = hist
        ! restraining potential of window i at coordinate xi_j
        w(i, :) = dexp(-beta * v)
    end do
#if defined(_OPENMP)
    !$omp end parallel do
#endif

    call unbias(w, p_biased)

end program

subroutine unbias(w, p_biased)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    implicit none

    integer i, j, k
    real*8, intent(in) :: w(nw, n), p_biased(nw, n)
    real*8 eps, fi_old, pmin, tmp1, tmp2, wtime, st, ft
    real*8 numerator, denominator, fi_new, invbeta
    real*8 p_unbiased(n), fi(nw), tmp(n)
    character(len=30) ::  date
    logical ex

    fi = 1.0 ! initial guess of unity for fi = dexp(-fi*beta)
    eps = tol
    k = 0
#if defined(_OPENMP)
    st = omp_get_wtime()
#else
    call cpu_time(st)
#endif
    do while(eps >= tol)
#if defined(_OPENMP)
        !$omp parallel do &
        !$omp private(j, i, denominator, numerator, tmp1) &
        !$omp shared(ni, w, fi, p_biased, p_unbiased)
#endif
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
            if(tmp1 == 0) tmp1 = 1.d-15
            p_unbiased(j) = tmp1
        end do
#if defined(_OPENMP)
        !$omp end parallel do
#endif

        !compute new fi based and the old use eps=sum{(1-fi_new/fi_old)**2}
        eps = 0.d0
#if defined(_OPENMP)
        !$omp parallel do &
        !$omp private(j, i, fi_new, fi_old) &
        !$omp shared(w, fi, p_unbiased) &
        !$omp reduction(+: eps)
#endif
        do i = 1, nw
            fi_old = fi(i)
            fi_new = 0.d0
            do j = 1, n
                fi_new = fi_new + w(i, j) * p_unbiased(j)
            end do
            fi(i) = fi_new
            eps = eps + (1.d0 - fi_new / fi_old) ** 2
        end do
#if defined(_OPENMP)
        !$omp end parallel do
#endif
        k = k + 1
        if(mod(k, 10000) == 1) then
            write(*, *) k, eps
            inquire(file="../exit", exist=ex)
            if(ex) call stopgm(nw, 'soft exit')
        end if
    end do
#if defined(_OPENMP)
    ft = omp_get_wtime()
#else
    call cpu_time(ft)
#endif
    wtime = ft - st
    !open(10,file='free_ener.dat', access='append')
    write(10, '(a,f6.2,a)') '# Unbiasing took ', wtime, 's'
    write(*, '(a, i9, a)') 'converged after ', k, ' loops'

    !find free energy shift
    invbeta = 1.d0 / beta
    pmin = invbeta * dlog(maxval(p_unbiased))
    if(symm) then
        k = maxloc(p_unbiased, dim=1)
        j = n / 2
        if(k > j) then
            tmp = [p_unbiased(n:n-j+1:-1), p_unbiased(j+1:n)]
        else
            tmp = [p_unbiased(1:j), p_unbiased(j+1:1:-1)]
        end if
    else
        tmp = p_unbiased
    end if

    !print pmf
    !open(10,file='free_ener.dat', access='append')
    do j = 1, n
        tmp1 = xbin(j) / fac
        tmp2 = (-invbeta * dlog(tmp(j)) + pmin) * fac1
        write(10, '(2f12.7)') tmp1, tmp2
    end do
    call movefile(nw, -1)

    call fdate(date)
    write(10, '(2a)') '# program ended on ', date
    close(10)
    write(*, *) 'pmf printed to free_ener.dat'
end subroutine
