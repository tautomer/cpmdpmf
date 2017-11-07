subroutine prob(i, hist, v)
    use global
    implicit none

    integer, intent(in) :: i
    real*8, intent(out) :: hist(n), v(n)
    integer j, k, m, nc, nat, trash, bin
    real*8 x
    real*8 tmp, tmp1, tmp2, tmp3, dist, dist1, dist2
    character(len=12) :: str_i
    character(len=20) :: input, output, test  ! for file reading purposesa
    character(len=*), parameter :: prefix = "traj_", test_sub = ".test", &
                                   & out_sub = ".output"

    write(str_i, *) i
    if(nb.eq.1) then
        input = 'CONSTRAINT'
    else
        input = 'TRAJECTORY'
    end if
    output = trim(prefix // adjustl(str_i)) // out_sub
    test = trim(prefix // adjustl(str_i)) // test_sub
    write(*, *) input, output
    open(unit=20,file=input,status='old')
    open(unit=21,file=output,status='unknown')
    open(unit=22,file=test,status='unknown')

    hist = 0
    if(nb.eq.1) then
        call read_traj(i, hist)
    else
        call read_rpmd_traj(i, hist)
    end if

    call get_biased(i, hist, v)

    return
end subroutine

subroutine read_traj(i, hist)
    use global
    implicit none

    integer k, bin, junk
    real*8 dist, tmp
    integer, intent(in) :: i
    real*8, intent(out) :: hist(n)

    do k = 1, ncut
        read(20, *)
    end do

    do k = ncut + 1, nsteps(i)
        read(20, *) junk, junk, tmp, dist
        dist = dist - xmin + tmp
        dist = dist / wbin
        if((dist.ge.n).or.(dist.lt.0))  cycle
        bin = 1 + dint(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        write(22,*) hist(bin), bin
    end do
    close(20)
    close(22)
    return
end subroutine

subroutine read_rpmd_traj(i, hist)
    use global
    implicit none

    integer j, k, l, m, bin, junk, nline
    real*8 r(natom, 3), d(3, 3), dist, tmp, l1, l2
    integer, intent(in) :: i
    real*8, intent(out) :: hist(n)

    nline = nb * natom
    do k = 1, ncut
        do j = 1, nline
            read(20, *)
        end do
    end do

    do j = ncut + 1, nsteps(i)
        d = 0
        do k = 1, nb
            do l = 1, natom
                read(20, *) junk, r(l, :)
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
        write(22,*) hist(bin), bin
    end do
    close(20)
    close(22)
    return
end subroutine

subroutine get_biased(i, hist, v)
    use global
    implicit none

    integer j
    real*8 tmp, k
    integer, intent(in) :: i
    real*8, intent(inout) :: hist(n)
    real*8, intent(out) :: v(n)

    k = ks(i) / 2
    !     compute normalized distribution and biasing potential
    write(21,'(a)') '# coordinate     potential     probability'
    do j = 1, n
        hist(j) = hist(j) / ni(i) ! normalized probality at tmp3
        tmp = xbin(i) - xi(i)
        v(j) = k * tmp ** 2  ! biasing window potential
        !   v(i)    = v(i)*fac1**2  ! convert potential to a.u.
        write(21, '(4f22.7)') xbin(j), v(j), hist(j)
    end do
    close(21)
    return
end subroutine
