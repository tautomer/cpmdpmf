subroutine run_wham(prim_pmf)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use wham
    implicit none

    integer :: i
    ! w(:,:)          biasing potential
    ! p_biased(:,:)   biased distribution
    ! p_unbiased(:)
    ! fi(:)           fi = exp(-fi*beta)
    ! hist(:)
    ! v(:)
    real*8, allocatable :: hist(:), v(:), w(:,:), p_biased(:,:)
    real*8, intent(out) :: prim_pmf(n)

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

    call unbias(w, p_biased, prim_pmf)

    return
end subroutine

subroutine prob(i, hist, v)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use wham
    implicit none

    integer :: uin
    character(len=11) :: flnm
    character(len=20) :: input, output, test
    integer, intent(in) :: i
    real*8, intent(out) :: hist(n), v(n)

    if(nb == 1) then
        flnm = "/CONSTRAINT"
    else 
        flnm = "/TRAJECTORY"
    end if
    write(input, "(3a)") "../", trim(dir(i)), flnm
    write(output, "(2a)") trim(dir(i)), ".output"
!    write(test, "(2a)") trim(dir(i)), ".test"
    write(*, *) i, output
#if defined(_OPENMP)
    uin = 20 + omp_get_thread_num()
#else
    uin = 20
#endif
    open(unit=uin, file=input, status="old")
    uin = uin + 2 * nw
    open(unit=uin, file=output, status="unknown")
    uin = uin - nw
    !open(unit=uin, file=test, status="unknown")

    hist = 0
    if(nb == 1) then
        call read_traj(i, hist, uin)
    else
        call read_rpmd_traj(i, hist, uin)
    end if
    call get_biased(i, hist, v, uin)

    return
end subroutine

subroutine read_traj(i, hist, udebug)
    use global
    use wham
    implicit none

    integer :: k, bin, junk, uin, ioerr
    real*8 :: dist
    integer, intent(in) :: i, udebug
    real*8, intent(out) :: hist(n)

    uin = udebug - nw
    do k = 1, ncut
        read(uin, *, iostat=ioerr)
        if(ioerr /= 0) call stopgm("No enough data in window ", dir(i))
    end do

    outer: do
        read(uin, *, iostat=ioerr) junk, junk, dist, dist
        if(ioerr > 0) then
            call stopgm("Error in reading data from window ", dir(i))
        else if(ioerr < 0) then
            exit
        end if
        dist = dist - xmin + xi(i)
        dist = dist / wbin
        if((dist > n).or.(dist < 0))  cycle
        bin = 1 + int(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        do bin = 1, nskip - 1
            read(uin, *, iostat=ioerr)
            if(ioerr < 0) exit outer
        end do
        !write(udebug, *) hist(bin), bin
    end do outer
    close(uin)
    !close(udebug)
    return
end subroutine

subroutine read_rpmd_traj(i, hist, udebug)
    use global
    use wham
    implicit none

    integer :: j, k, bin, junk, nloop, uin, counter, ioerr
    real*8 :: coor(natom, 3), proj, diff, dist
    integer, intent(in) :: i, udebug
    real*8, intent(out) :: hist(n)

    uin = udebug - nw
    nloop = nb * natom * ncut
    do j = 1, nloop
        read(uin, *, iostat=ioerr)
        if(ioerr /= 0) call stopgm("No enough data in window ", dir(i))
    end do

    nloop = nb * natom * (nskip - 1)
    outer: do
        dist = 0
        do j = 1, nb
            do k = 1, natom
                read(uin, *, iostat=ioerr) junk, coor(k, :)
                if(ioerr > 0) then
                    call stopgm("Error in reading data from window ", dir(i))
                else if(ioerr < 0) then
                    exit outer
                end if
            end do
            if (tproj) then
                dist = dist + proj(coor(ind(:), :))
            else 
                dist = dist + diff(coor(ind(:), :))
            end if
        end do
        dist = dist / nb
        dist = dist - xmin
        dist = dist / wbin
        if((dist < 0).or.(dist > n))  cycle
        bin = 1 + int(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        !write(udebug,*) hist(bin), bin
        do j = 1, nloop
            read(uin, *, iostat=ioerr)
            if(ioerr < 0) exit outer
        end do
        if(nread > ncut) then
            counter = counter + nskip
            if(counter >= nread) exit
        end if
    end do outer
    close(uin)
    close(udebug)

    return
end subroutine

subroutine get_biased(i, hist, v, udebug)
    use global
    use wham
    implicit none

    integer :: j, uout
    real*8 :: tmp, k
    integer, intent(in) :: i, udebug
    real*8, intent(inout) :: hist(n)
    real*8, intent(out) :: v(n)

    uout = udebug + nw
    k = ks(i) / 2
    ! compute normalized distribution and biasing potential
    write(uout,"(a)") "# coordinate     probability     potential"
    do j = 1, n
        hist(j) = hist(j) / ni(i) ! normalized probality at tmp3
        tmp = xbin(j) - xi(i)
        v(j) = k * tmp ** 2  ! biasing window potential
        write(uout, "(3f20.7)") xbin(j), hist(j), v(j)
    end do
    close(uout)
    return
end subroutine

subroutine unbias(w, p_biased, prim_pmf)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use wham
    implicit none

    integer :: i, j, k
    real*8 :: eps, fi_old, tmp1, wtime, st, ft
    real*8 :: numerator, denominator, fi_new, pmin
    real*8 :: p_unbiased(n), fi(nw) 
    logical :: ex
    real*8, intent(in) :: w(nw, n), p_biased(nw, n)
    real*8, intent(out) :: prim_pmf(n)

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
            if(ex) call stopgm("Soft exit")
        end if
    end do
#if defined(_OPENMP)
    ft = omp_get_wtime()
#else
    call cpu_time(ft)
#endif
    wtime = ft - st
    !open(10,file="free_ener.dat", access="append")
    write(10, "(a,f6.2,a)") "# Unbiasing took ", wtime, "s"
    write(*, "(a, i9, a)") "converged after ", k, " loops"

    prim_pmf = -kt * dlog(p_unbiased)
    pmin = minval(prim_pmf)
    prim_pmf = prim_pmf - pmin

    return
end subroutine
