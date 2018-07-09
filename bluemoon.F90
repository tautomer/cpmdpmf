subroutine run_blue_moon(prim_pmf)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use bluemoon, only: dx
    implicit none

    integer :: i
    character(len=30) :: msg
    real*8 :: dadx(n)
    real*8, intent(out) :: prim_pmf(n)

#if defined(_OPENMP)
    !$omp parallel do &
    !$omp private(i) &
    !$omp shared(dadx)
#endif
    do i = 1, nw ! loop over windows
        call get_da_dx(i, dadx(i), msg)
        if(msg /= " ") call stopgm(msg)
    end do
#if defined(_OPENMP)
    !$omp end parallel do

    prim_pmf(1) = 0.0
    do i = 2, nw
        prim_pmf(i) = prim_pmf(i-1) + (dadx(i) + dadx(i-1)) * dx(i-1) * fac1
    end do
#endif

    return
end subroutine

subroutine get_da_dx(i, dadx, msg)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use bluemoon
    implicit none

    integer :: ucons, utraj
    character(len=20) :: fcons, ftraj
    character(len=*), parameter :: cons = "/CONSTRAINT", traj = "/TRAJECTORY"
    integer, intent(in) :: i
    real*8, intent(out) :: dadx
    character(len=*), intent(out) :: msg

    write(fcons, '(3a)') "../", trim(dir(i)), cons
    write(ftraj, '(3a)') "../", trim(dir(i)), traj
#if defined(_OPENMP)
    ucons = 20 + omp_get_thread_num()
#else
    utraj = 20
#endif
    open(unit=ucons, file=fcons, status='old')
    utraj = ucons + nw
    open(unit=utraj, file=ftraj, status='unknown')

    write(*, *) i
    call read_data(i, ucons, dadx, msg)

    return
end subroutine

subroutine read_data(i, ucons, dadx, msg)
    use global
    use bluemoon
    implicit none

    integer :: j, k, junk, utraj, err, nloop, counter, tmp
    real*8 :: lambda, coor(natom, 3), g, sum_denom, sum_numer, sum_lambda
    real*8 :: r(3, 3), dsdx(3, 3), d2sdx2(3, 3, 3, 3), invzr, test
    integer, intent(in) :: i, ucons
    real*8, intent(out) :: dadx
    character(len=*), intent(out) :: msg

    msg = " "
    utraj = ucons + nw
    nloop = nb * ncut
    do j = 1, nloop
        read(ucons, *, iostat=err)
        if(err /= 0) call stopgm('No enough data in window ', dir(i))
        do k = 1, natom
            read(utraj, *, iostat=err)
            if(err /= 0) call stopgm('No enough data in window ', dir(i))
        end do
    end do

    nloop = nb * (nskip - 1)
    sum_numer = 0
    sum_denom = 0
    sum_lambda = 0
    counter = 0
    tmp = 0
    outer: do
        r = 0
        do j = 1, nb
            read(ucons, *, iostat=err) junk, junk, lambda
            if(err > 0) then
                call stopgm('Error in reading data from window ', dir(i))
            else if(err < 0) then
                exit outer
            end if
            do k = 1, natom
                read(utraj, *, iostat=err) junk, coor(k, :)
                if(err > 0) then
                    call stopgm('Error in reading data from window ', dir(i))
                else if(err < 0) then
                    exit outer
                end if
            end do
            r(1, :) = r(1, :) + coor(ind(1), :)
            r(2, :) = r(2, :) + coor(ind(2), :)
            r(3, :) = r(3, :) + coor(ind(3), :)
        end do
        !sum_lambda = sum_lambda + lambda
        !tmp = tmp + 1
        r = r / nb
        call get_zr(r, dsdx, invzr)
        call get_2nd_deriv(r, d2sdx2)
        call get_g(dsdx, d2sdx2, invzr, g)
        call sum_numer_denom(invzr, lambda, g, sum_numer, sum_denom)
        !write(50+i, *) tmp, sum_denom / tmp
        !write(100+i, *) tmp, sum_lambda / tmp
        do j = 1, nloop
            read(ucons, *, iostat=err)
            if(err < 0) exit outer
            do k = 1, natom
                read(utraj, *, iostat=err)
                if(err < 0) exit outer
            end do
        end do
        if(nread > ncut) then
            counter = counter + nskip
            if(counter >= nread) exit
        end if
    end do outer
    close(ucons)
    close(utraj)
    dadx = sum_numer / sum_denom !* facda

    return
end subroutine

subroutine get_zr(r, dsdx, invzr)
    use global
    use bluemoon
    implicit none

    integer :: i, j
    real*8 :: newr(3, 3), proj, diff, f1, f2, zr
    real*8, intent(in) :: r(3, 3)
    real*8, intent(out) :: dsdx(3, 3), invzr

    newr = r
    if(tdiff) then
        do i = 1, 3
            do j = 1, 3        
                newr(i, j) = r(i, j) + h
                f1 = diff(newr)
                newr(i, j) = r(i, j) - h
                f2 = diff(newr)
                dsdx(i, j) = (f1 - f2) * h2
                newr(i, j) = r(i, j)
            end do      
        end do
    else
        do i = 1, 3
            do j = 1, 3        
                newr(i, j) = r(i, j) + h
                f1 = proj(newr)
                newr(i, j) = r(i, j) - h
                f2 = proj(newr)
                dsdx(i, j) = (f1 - f2) * h2
                newr(i, j) = r(i, j)
            end do      
        end do
    end if
    zr = 0
    do i = 1, 3
        do j = 1, 3
            zr = zr + invm(i) * dsdx(i, j) ** 2
        end do
    end do
    invzr = 1 / zr

    return
end subroutine

subroutine get_2nd_deriv(r, d2sdx2)
    use global, only: tdiff
    use bluemoon
    implicit none

    integer :: i, j, k, l, m, n
    real*8 :: newr(3, 3), proj, diff, f(4)
    real*8, intent(in) :: r(3, 3)
    real*8, intent(out) :: d2sdx2(3, 3, 3, 3)

    newr = r
    if(tdiff) then
        do i = 1, 3
            do k = 1, 3
                do j = 1, 3
                    do l = 1, 3
                        do m = 1, 4
                            n = m / 2 + 1
                            newr = r
                            newr(i, k) = newr(i, k) - (-1) ** m * h
                            newr(j, l) = newr(j, l) - (-1) ** n * h
                            f(m) = diff(newr)
                            newr(i, k) = newr(i, k)
                            newr(j, l) = newr(j, l)
                        end do
                        d2sdx2(i, k, j, l) = (f(1) + f(2) - f(3) - f(4)) * h42
                    end do
                end do
            end do
        end do
    else
        do i = 1, 3
            do k = 1, 3
                do j = 1, 3
                    do l = 1, 3
                        do m = 1, 4
                            n = m / 2 + 1
                            newr = r
                            newr(i, k) = newr(i, k) - (-1) ** m * h
                            newr(j, l) = newr(j, l) - (-1) ** n * h
                            f(m) = proj(newr)
                            newr(i, k) = newr(i, k)
                            newr(j, l) = newr(j, l)
                        end do
                        d2sdx2(i, k, j, l) = (f(1) + f(2) - f(3) - f(4)) * h42
                    end do
                end do
            end do
        end do
    end if

    return
end subroutine

subroutine get_g(dsdx, d2sdx2, invzr, g)
    use bluemoon
    implicit none

    integer i, j, k, l
    real*8, intent(in) :: invzr, dsdx(3, 3), d2sdx2(3, 3, 3, 3)
    real*8, intent(out) :: g

    g = 0
    do i = 1, 3
        do k = 1, 3
            do j = 1, 3
                do l = 1, 3
                    g = g + invmij(i, j) * dsdx(i, k) * d2sdx2(i, k, j, l) * &
                        dsdx(j, l)
                end do
            end do
        end do
    end do
    g = g * invzr ** 2
    return
end subroutine


subroutine sum_numer_denom(invzr, lambda, g, sum_numer, sum_denom)
    use global, only: kt
    use bluemoon
    implicit none

    real*8 sqrtzr
    real*8, intent(in) :: invzr, lambda, g
    real*8, intent(inout) :: sum_numer, sum_denom

    sqrtzr = dsqrt(invzr)
    sum_numer = sum_numer + sqrtzr * (lambda + kt * g)
    sum_denom = sum_denom + sqrtzr

    return
end subroutine