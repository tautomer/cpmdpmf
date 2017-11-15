subroutine read_conf()
    use global
    implicit none

    logical ex

    open(unit=20,file="config", action="read")
    read(20, *) nw
    allocate(dir(nw), xi(nw), nsteps(nw), ni(nw), ks(nw))
    read(20, *) wbin ! unit in angstrom
    read(20, *) ncut
    read(20, *) temp
    read(20, *) tol
    read(20, *) dir
    close(20)
    !!! inquire works differently in ifort and gfortran
    !!! thinking to use a $FC variable in makefile to address this
    inquire(directory="./debug", exist=ex)
    if (ex) then
        write(*, *) "Dir debug already exists. Cleaning it up."
    !!! if previous job wasn't terminated finished, doing so may lose files
    !!! currently I can back everything up once
    !!! might be solved by inquire
        call system("rm -f debug/*")
    else
        call system("mkdir -p debug")
    end if

end subroutine

subroutine folderloop(sgn)
    use global
    implicit none

    integer i
    integer, intent(in) :: sgn
    character(len=12) :: path, rmtmp
    character(len=70) :: getcfg
    character(len=200) :: getmol

    getcfg = "grep DIF inp-2 | cut -d' ' -f4-8 > tmpin; sed -n '/MAXS/{n;p;}'&
             & inp-2 >> tmpin"
    getmol = "m=$(grep -c TROT inp-2); if [[ $m -eq 0 ]]; then nb=1; else nb=$(&
             &sed -n '/TROT/{n;p;}' inp-2); fi; nat=$(wc -l $(ls GEOMETRY*.xyz &
             &| head -1) | cut -d' ' -f1); echo $nb $nat >> tmpin"
    rmtmp = "rm -f tmpin"
    call getcwd(rootdir)
    open(unit=11, file='xi')
    read(11, *) xi
    close(11)
    ks(:51) = 0.1
    ks(52:) = 0.3
    nb = 1
    natom = 1
    nsteps = 3d4
    do i = 1, nw
        path = trim("./" // adjustl(dir(i)))
        call chdir(path)
!        if(sgn.eq.1) then
!            call system(getcfg)
!            if(i.eq.1) call system(getmol)
!            open(unit=11, file='tmpin')
!            read(11, *) ind, xi(i), ks(i)
!            read(11, *) nsteps(i)
!            if(i.eq.1) read(11, *) nb, natom
!            close(11)
!        end if
        call movefile(i, sgn)
        call chdir(rootdir)
    end do
end subroutine

subroutine init_param(date)
    use omp_lib
    use global
    implicit none
    integer i, ind1, ind2
    real*8 x(nw)
    character(len=30), intent(in) ::  date

    wbin = wbin * 1
    ! in case the xi array is not sorted
    ! adding a sort function?
    ind1 = maxloc(xi, 1)
    x = [xi(:ind1-1), xi(ind1+1:)]
    ind2 = maxloc(x, 1)
    xmax = (3 * xi(ind1) - x(ind2)) * 1 * 0.5
    ind1 = minloc(xi, 1)
    x = [xi(:ind1-1), xi(ind1+1:)]
    ind2 = minloc(x, 1)
    xmin = (3 * xi(ind1) - x(ind2)) * 1 * 0.5
    n = dint((xmax - xmin) / wbin) ! total number of bins
    allocate(xbin(n))
    beta = 1.d0 / kb / temp
    ni = 0
    do i = 1, n
        xbin(i) = xmin + wbin * (i - 0.5)
    end do
    !!! need a folder to store data of each window
    open(10, file='free_ener.dat', status='unknown')
    write(10, *) '# number of beads: ', nb
    write(10, *) '# beta: ', beta
    write(10, *) '# number of windows: ', nw
    write(10, *) '# program started on ', date
    write(10, '(a,i8)' ) &
        ' # the number of processors available: ', omp_get_num_procs ()
    write(10, '(a,i8)' ) &
        ' # the number of threads available: ', omp_get_max_threads ()
    call chdir('./debug')
end subroutine

subroutine movefile(i, sgn)
    use global

    integer, intent(in) :: i, sgn
    character(len=12) :: str_i, input
    character(len=30) :: output
    character(len=70) :: cmd

        if(nb.eq.1) then
            input = ' CONSTRAINT'
        else
            input = ' TRAJECTORY'
        end if
        write(str_i, *) i
        output = trim(' ../debug/traj_' // adjustl(str_i)) // '.dat'
        if(sgn.eq.1) then
            cmd = trim('mv ' // adjustl(input)) // output
        else
            cmd = trim('mv ' // adjustl(output)) // input
        end if
        call system(cmd)
end subroutine
