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

    ! the output normalized probability distribution
    !$OMP PARALLEL DO
    do j = 1, n               ! initial histogram values to zero
        hist(j) = 0
    end do
    !$OMP END PARALLEL DO

    !nc = 0 ! initialize total number of snapshots to zero

    !allocate(pos(3,3))

    !     compute normalized distribution and biasing potential
    write(21,'(a)') '# coordinate     potential     probability'
    do j = 1, n
        hist(j) = hist(j) / ni(i) ! normalized probality at tmp3
        tmp2 = xbin(i) - xi(i)
        v(j) = ks / 2 * tmp2 ** 2  ! biasing window potential
        !   v(i)    = v(i)*fac1**2  ! convert potential to a.u.
        write(21,'(4f22.7)') tmp3, v(j), hist(j)
    end do
    close(21)
    return
end subroutine

subroutine read_traj(i, hist)
    implicit none

    integer k, bin, junk
    real*8 dist, tmp
    integer, intent(in) :: i
    real*8, intent(out) :: hist(n)

    do k = 1, ncut
        read(20, *)
    end do

    do k = ncut + 1, nsteps(i)
        read(20, *) j, junk, tmp, dist
        dist = dist - xmin + tmp
        dist = dist / wbin
        if((dist.ge.n).or.(dist.lt.0))  cycle
        bin = 1 + dint(dist) ! locate bin
        ni(i) = ni(i) + 1
        hist(bin) = hist(bin) + 1 ! place in appropriate bin
        write(22,*) hist(bin),bin
    end do
    close(20)
    close(22)
