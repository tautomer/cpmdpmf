subroutine prob(i, hist, v)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    implicit none

    integer, intent(in) :: i
    real*8, intent(out) :: hist(n), v(n)
    integer uin
    character(len=12) :: str_i
    character(len=20) :: input, output, test
    character(len=*), parameter :: prefix = "traj_", test_sub = ".test", &
                                   & out_sub = ".output"

    write(str_i, *) i
    input = trim(prefix // adjustl(str_i)) // '.dat'
    output = trim(prefix // adjustl(str_i)) // out_sub
    test = trim(prefix // adjustl(str_i)) // test_sub
    write(*, *) i, output
#if defined(_OPENMP)
    uin = 20 + omp_get_thread_num()
#else
    uin = 20
#endif
    open(unit=uin, file=input, status='old')
    uin = uin + 2 * nw
    open(unit=uin, file=output, status='unknown')
    uin = uin - nw
    !open(unit=uin, file=test, status='unknown')

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
    implicit none

    integer k, bin, junk, uin, err
    real*8 dist
    integer, intent(in) :: i, udebug
    real*8, intent(out) :: hist(n)

    uin = udebug - nw
    do k = 1, ncut
        read(uin, *, iostat=err)
        if(err /= 0) call stopgm(nw, 'no enough data in window ', dir(i))
    end do

    outer: do
        read(uin, *, iostat=err) junk, junk, dist, dist
        if(err > 0) then
            call stopgm(nw, 'error in reading data from window ', dir(i))
        else if(err < 0) then
            exit
        end if
        dist = dist - xmin + xi(i)
        dist = dist / wbin
        if((dist > n).or.(dist < 0))  cycle
        bin = 1 + dint(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        do bin = 1, nskip - 1
            read(uin, *, iostat=err)
            if(err < 0) exit outer
        end do
        !write(udebug, *) hist(bin), bin
    end do outer
    close(uin)
    close(udebug)
    return
end subroutine

subroutine read_rpmd_traj(i, hist, udebug)
    use global
    implicit none

    integer j, k, bin, junk, nline, uin, err
    real*8 r(natom, 3), xy(3), xz(3), dxz, r2, d
    real*8 dist, l1, l2
    integer, intent(in) :: i, udebug
    real*8, intent(out) :: hist(n)

    uin = udebug - nw
    nline = nb * natom * ncut
    do j = 1, nline
        read(uin, *, iostat=err)
        if(err /= 0) call stopgm(nw, 'no enough data in window ', dir(i))
    end do

    nline = nb * natom * (nskip - 1)
    outer: do
        dist = 0
        do j = 1, nb
            do k = 1, natom
                read(uin, *, iostat=err) junk, r(k, :)
                if(err > 0) then
                    call stopgm(nw, 'error in reading data from window ', &
                                dir(i))
                else if(err < 0) then
                    exit outer
                end if
            end do
            xy = r(ind(2), :) - r(ind(1), :)
            xz = r(ind(3), :) - r(ind(1), :)
            r2 = 0
            r2 = r2 + xz(1) ** 2
            r2 = r2 + xz(2) ** 2
            r2 = r2 + xz(3) ** 2
            dxz = dsqrt(r2)
            d = 0
            d = d + xy(1) * xz(1)
            d = d + xy(2) * xz(2)
            d = d + xy(3) * xz(3)
            dist = dist + d / dxz
        end do
        dist = dist / nb
        dist = dist - xmin
        dist = dist / wbin
        if((dist < 0).or.(dist > n))  cycle
        bin = 1 + dint(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        !write(udebug,*) hist(bin), bin
        do j = 1, nline
            read(uin, *, iostat=err)
            if(err < 0) exit outer
        end do
    end do outer
    close(uin)
    close(udebug)
    return
end subroutine

subroutine get_biased(i, hist, v, udebug)
    use global
    implicit none

    integer j, uout
    real*8 tmp, k
    integer, intent(in) :: i, udebug
    real*8, intent(inout) :: hist(n)
    real*8, intent(out) :: v(n)

    uout = udebug + nw
    k = ks(i) / 2
    ! compute normalized distribution and biasing potential
    write(uout,'(a)') '# coordinate     probability     potential'
    do j = 1, n
        hist(j) = hist(j) / ni(i) ! normalized probality at tmp3
        tmp = xbin(j) - xi(i)
        v(j) = k * tmp ** 2  ! biasing window potential
        write(uout, '(3f20.7)') xbin(j), hist(j), v(j)
    end do
    close(uout)
    return
end subroutine
