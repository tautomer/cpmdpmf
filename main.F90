program recons_pmf ! change this name
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    implicit none

    real*8, allocatable :: prim_pmf(:)
    call fdate(date)
    call read_conf()
    call folder_loop()
    call init_param()
    allocate(prim_pmf(n))

    if (twham) then
        call run_wham(prim_pmf)
    else
        call run_blue_moon(prim_pmf)
    end if

    call write_pmf(prim_pmf)
end program

subroutine write_pmf(prim_pmf)
    use global
    use bluemoon
    use wham
    implicit none

    integer :: j, k
    real*8 :: tmp(n), tmp1
    real*8, intent(in) :: prim_pmf(n)

    if(symm) then
        k = minloc(prim_pmf, dim=1)
        j = n / 2
        if(k > j) then
            tmp = [prim_pmf(n:n-j+1:-1), prim_pmf(j+1:n)]
        else
            tmp = [prim_pmf(1:j), prim_pmf(j+1:1:-1)]
        end if
    else
        tmp = prim_pmf
    end if

    !print pmf
    !open(10,file='free_ener.dat', access='append')
    do j = 1, n
        tmp1 = xbin(j) / fac
        write(10, '(2f12.7)') tmp1, tmp(j)
    end do

    call fdate(date)
    write(10, '(2a)') '# program ended on ', date
    close(10)
    write(*, *) 'pmf printed to free_ener.dat'

end subroutine
