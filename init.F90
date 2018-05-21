subroutine read_conf()
    use global
    implicit none

    logical ex

    open(unit=20, file="config", action="read")
    read(20, *) nw
    allocate(dir(nw), xi(nw), nsteps(nw), ni(nw), ks(nw))
    read(20, *) wbin ! unit in angstrom
    read(20, *) ncut
    read(20, *) nskip
    read(20, *) temp
    read(20, *) tol
    read(20, *) symm
    read(20, *) dir
    close(20)
!10  call stopgm(0, 'config file error')
    !!! inquire works differently in ifort and gfortran
    !!! thinking to use a $FC variable in makefile to address this
    inquire(directory="./debug", exist=ex)
    if (ex) then
        write(*, *) "Dir debug already exists. Cleaning it up."
    !!! if previous job wasn't terminated finished, doing so may lose files
    !!! currently I can back everything up once
    !!! might be solved by inquire
        call system("rm -f debug/* exit")
    else
        call system("mkdir -p debug")
    end if

end subroutine

subroutine folderloop()
    use global
    implicit none

    integer i
    character(len=12) :: path

    call getcwd(rootdir)
    do i = 1, nw
        path = trim("./" // adjustl(dir(i)))
        call chdir(path)
        call read_tmp(i)
        call chdir(rootdir)
    end do
end subroutine

subroutine read_tmp(i)
    use ifport
    use global

    integer, intent(in) :: i
    integer msg

    character(len=80) :: getcfg
    character(len=250) :: getmol

    getcfg = "grep PROJ inp-2 | cut -d' ' -f4-8 > tmpin; sed -n '/MAXS/{n;p;}'&
             & inp-2 >> tmpin"
    getmol = "m=$(grep -c INTEG inp-2); if [[ $m -eq 0 ]]; then nb=1; else nb=&
             $(sed -n '/TROT/{n;p;}' inp-2); fi; geo=$(ls GEOMETRY*.xyz | hea&
             d -1); [[ -z $geo ]] && exit 1; nat=$(wc -l $geo|cut -d' ' -f1); &
             echo $nb $(($nat-2)) >> tmpin"
    msg = system(getcfg)
    if(msg /= 0) call stopgm(0, 'inp-2 not found in window ', dir(i))
    if(i == 1) then
        msg = system(getmol)
        if(msg /= 0) call stopgm(0, 'GEOMETRY*.xyz not found in window ', &
                                 dir(i))
    end if
    open(unit=11, file='tmpin')
    read(11, *) ind, xi(i), ks(i)
    read(11, *) nsteps(i)
    if(i == 1) read(11, *) nb, natom
!20  call stopgm(0, 'reading inp-2 error')
    close(11, status='delete')
end subroutine

subroutine init_param(date)
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    implicit none
    integer i, ind1, ind2
    real*8 x(nw)
    character(len=30), intent(in) ::  date

    wbin = wbin * fac
    xi = xi * fac
    !nsteps = (nsteps - ncut) / nskip
    ! in case the xi array is not sorted
    ! adding a sort function?
    ind1 = maxloc(xi, 1)
    x = [xi(:ind1-1), xi(ind1+1:)]
    ind2 = maxloc(x, 1)
    xmax = (3 * xi(ind1) - x(ind2)) * 0.5
    ind1 = minloc(xi, 1)
    x = [xi(:ind1-1), xi(ind1+1:)]
    ind2 = minloc(x, 1)
    xmin = (3 * xi(ind1) - x(ind2)) * 0.5
    n = dint((xmax - xmin) / wbin) ! total number of bins
    allocate(xbin(n))
    beta = 1.d0 / kb / temp
    ni = 0
#if defined(_OPENMP)
    !$omp parallel do &
    !$omp private(i) &
    !$omp shared(xbin, xmin, wbin)
#endif
    do i = 1, n
        xbin(i) = xmin + wbin * (i - 0.5)
    end do
#if defined(_OPENMP)
    !$omp end parallel do
#endif
    open(10, file='free_ener.dat', status='unknown')
    write(10, '(a,i8)') '# number of beads: ', nb
    write(10, '(a,f8.2)') '# beta: ', beta
    write(10, '(a,i8)') '# number of windows: ', nw
    write(10, '(2a)') '# program started on ', date
#if defined(_OPENMP)
    write(10, '(a,i8)') &
        '# the number of processors available: ', omp_get_num_procs ()
    write(10, '(a,i8)') &
        '# the number of threads available: ', omp_get_max_threads ()
#else
    write(10, '(a)' ) '# OMP disabled'
#endif
    call chdir('./debug')
end subroutine

recursive subroutine movefile(n, sgn)
    use ifport
    use global, only : dir, nb

    integer, intent(in) :: n, sgn
    integer i, err
    character(len=12) :: path, idx
    character(len=30) :: inpnm, input, output
    character(len=70) :: cmd

    if(nb == 1) then
        inpnm = '/CONSTRAINT'
    else
        inpnm = '/TRAJECTORY'
    end if
    ! not implementing OMP
    ! error handling is depentent on sequence
    do i = 1, n
        write(path, *) dir(i)
        input = trim(' ../' // adjustl(path)) // inpnm
        write(idx, *) i
        output = trim(' traj_' // adjustl(idx)) // '.dat'
        if(sgn == 1) then
            cmd = trim('mv ' // adjustl(input)) // output
        else
            cmd = trim('mv ' // adjustl(output)) // input
        end if
        err = system(cmd)
        if(err /= 0) call stopgm(i-1, 'missing data from window ', dir(i))
    end do
end subroutine
