subroutine prob(i, hist, v)
    use omp_lib
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
    uin = 20 + omp_get_thread_num()
    open(unit=uin, file=input, status='old')
    uin = uin + 2 * nw
    open(unit=uin, file=output, status='unknown')
    uin = uin - nw
    !open(unit=uin, file=test, status='unknown')

    hist = 0
    if(nb.eq.1) then
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

    integer k, bin, junk, uin
    real*8 dist
    integer, intent(in) :: i, udebug
    real*8, intent(out) :: hist(n)

    uin = udebug - nw
    do k = 1, ncut
        read(uin, *)
    end do

    do k = ncut + 1, nsteps(i)
        read(uin, *) junk, dist
        dist = dist - xmin
        dist = dist / wbin
        if((dist.ge.n).or.(dist.lt.0))  cycle
        bin = 1 + dint(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        !write(udebug, *) hist(bin), bin
    end do
    close(uin)
    close(udebug)
    return
end subroutine

subroutine read_rpmd_traj(i, hist, udebug)
    use global
    implicit none

    integer j, k, l, m, bin, junk, nline, uin
    real*8 r(natom, 3), d(3, 3), dist, l1, l2
    integer, intent(in) :: i, udebug
    real*8, intent(out) :: hist(n)

    uin = udebug - nw
    nline = nb * natom
    do k = 1, ncut
        do j = 1, nline
            read(uin, *)
        end do
    end do

    do j = ncut + 1, nsteps(i)
        d = 0
        do k = 1, nb
            do l = 1, natom
                read(uin, *) junk, r(l, :)
            end do
            do l = 1, 3
                do m = 1, 3
                    d(l, m) = d(l, m) + r(ind(l), m)
                end do
            end do
        end do
        d = d / nb
        l1 = 0
        l2 = 0
        do l = 1, 3
            l1 = l1 + (d(1, l) - d(2, l)) ** 2
            l2 = l2 + (d(3, l) - d(2, l)) ** 2
        end do
        dist = sqrt(l1) - sqrt(l2)
        dist = dist - xmin
        dist = dist / wbin
        if((dist.ge.n).or.(dist.lt.0))  cycle
        bin = 1 + dint(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        !write(udebug,*) hist(bin), bin
    end do
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
    !     compute normalized distribution and biasing potential
    write(uout,'(a)') '# coordinate     potential     probability'
    do j = 1, n
        hist(j) = hist(j) / ni(i) ! normalized probality at tmp3
        tmp = xbin(j) - xi(i)
        v(j) = k * tmp ** 2  ! biasing window potential
        write(uout, '(4f20.7)') xbin(j), v(j), hist(j)
    end do
    close(uout)
    return
end subroutine
