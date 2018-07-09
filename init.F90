subroutine read_conf()
    implicit none

    logical :: ex
    character(len=30) :: config

    call get_command_argument(1, config)
    open(unit=20, file=config, action="read")
    call parse_conf()

#if defined(__INTEL_COMPILER)
    inquire(directory="./debug", exist=ex)
#else
    inquire(file="./debug", exist=ex)
#endif
    if (ex) then
        write(*, *) "Dir debug already exists. Cleaning it up."
        call system("rm -f debug/*")
    else
        call system("mkdir -p debug")
    end if
    call system("rm -f exit")

end subroutine

subroutine parse_conf()
    implicit none

    integer :: line = 0, ioerr = 0, flag = 0, pos
    character(len=100) :: buffer, label

    call def_val()
    do while(ioerr == 0)
        read(20, '(A)', iostat=ioerr) buffer
        if(ioerr == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1: pos)
            buffer = buffer(pos+1: )
            call read_buffer(label, buffer, line, flag)
        end if
    end do
    close(20)
    call fix_unspecd(flag)

end subroutine

subroutine def_val()
    use global
    use bluemoon
    use wham
    implicit none

    tbl = .false.
    twham = .false. 
    tdiff = .false. 
    tproj = .false. 
    nw = 0
    temp = 0
    tol = 0
    invm = 0
    ncut = -1
    nskip = -1
    nread = -1
    wbin = 0
end subroutine

subroutine read_buffer(label, buffer, line, flag)
    use global
    use bluemoon
    use wham
    implicit none

    integer :: ioerr
    integer, intent(in) :: line
    integer, intent(inout) :: flag
    character(len=100), intent(in) :: buffer, label
    select case (label)
    ! general settings
    ! need fix twham & tbl & tdiff & tproj
    case ("bluemoon")
        tbl = .true.
        twham = .false.
    case ("wham")
        tbl = .true.
        twham = .false.
    case ("differ")
        tdiff = .true.
        tproj = .false.
    case ("window")
        read(buffer, "(i0)", iostat=ioerr) nw
        if (ioerr > 0) then
            call stopgm("Number of windows should be an integer")
        else if (ioerr < 0) then
            call stopgm("Number of windows not sepcified")
        end if
        allocate(dir(nw), xi(nw), ni(nw), ks(nw))
    case ("temperature")
        read(buffer, *, iostat=ioerr) temp
        if (ioerr > 0) then
            call stopgm("Temperature should be a real number")
        else if (ioerr < 0) then
            write(*, "(a, f5.1, a)") "Temperature is specified, but no value &
            &provided. Use defualt value of ", 300, " k."
            temp = 300
        end if
    case ("directory")
        if (nw == 0) call stopgm("No of windows needed to be specified first.")
        read(buffer, *, iostat=ioerr) dir
        if (ioerr < 0) then
            call stopgm("No of dirs provided smaller than no of windows.")
        end if
    case ("symmetry")
        read(buffer, *, iostat=ioerr) symm
        if (ioerr > 0) then
            call stopgm("Symmetry of the system should be a logical variable.")
        else if (ioerr == 0) then
            flag = 1
        end if
    case ("step")
        read(buffer, "(3i0)", iostat=ioerr) ncut, nskip, nread
        if (ioerr > 0) then
            call stopgm("Step information should be an integer")
        else if (ioerr < 0) then
            write(*, "(a)") "Step is specified, but not all values are &
            &provided. Going to read all data."
            ncut = 0
            nskip = 0
            nread = 0
        end if
        if (nread > 0 .and. nskip >= nread) write(*, "(a)") "Only one step &
        &will be read. Double check values provided?"
    ! for WHAM
    case ("bin")
        read(buffer, *, iostat=ioerr) wbin
        if (ioerr > 0) then
            call stopgm("Bin size should be a real number.")
        else if (ioerr < 0) then
            write(*, "(a, f4.2, a)") "Bin is specified, but size is not. &
            &Use defualt value of ", 1d-2, " a.u."
            wbin = 1d-2
        end if
    case ("tolerance")
        read(buffer, *, iostat=ioerr) tol
        if (ioerr > 0) then
            call stopgm("SCF tolerance should be a real number.")
        else if (ioerr < 0) then
            write(*, "(a)") "Tolerance is specified, but no value &
            &provided. Use defualt value of 1d-9."
            tol = 1d-9
        end if
    ! for bluemoon
    case ("masses")
        read(buffer, *, iostat=ioerr) invm
        if (ioerr > 0) then
            call stopgm("Masses should be real numbers.")
        else if (ioerr < 0) then
            call stopgm("Not all masses are provided.")
        end if
    case default
        write(*, "(a, i0)") "Skipping invalid label at line ", line
    end select
    return

end subroutine

subroutine fix_unspecd(flag)
    use global
    use bluemoon
    use wham
    implicit none

    integer, intent(in) :: flag

    if (.not.twham .and. .not.tbl) call stopgm("Neither WHAM nor blue moon is &
    &specified.")
    if (temp == 0) then
        write(*, "(a)") "Temperature not specified. Use defualt 300K instead."
        temp = 300
    end if
    if (tol == 0) then
        write(*, "(a)") "Tolerance not specified. Use defualt 1d-9 instead."
        tol = 1d-9
    end if
    if (tbl .and. invm(1) == 0) call stopgm("Masses needed for blue moon.")
    if (ncut * nskip * nread < 0) then
        write(*, "(a)") "No information specified on steps. Going to read all &
        &data."
        ncut = 0
        nskip = 0
        nread = 0
    end if
    if (wbin == 0) then
        write(*, "(a)") "Bin size not specified. Use defualt 1d-2 a.u. instead."
        wbin = 1d-2
    end if
    if (flag == 1) then
        write(*, "(a)") "Symmetry not specified. The system is assumed to be &
        &asymmetric."
        symm = .false.
    end if

end subroutine

subroutine folder_loop()
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
#if defined(__INTEL_COMPILER)
    use ifport
#endif
    use global
    use wham
    implicit none

    integer, intent(in) :: i
    integer :: msg, ioerr

    character(len=250) :: getcfg*80, getmol

    ! the name of the input file 'inp-2' is hardcoded
    ! consider a better way to get all the information here
    getcfg = "grep DIFF inp-2 | cut -d' ' -f4-8 > tmpin; sed -n '/MAXS/{n;p;}'&
             & inp-2 >> tmpin"
    getmol = "m=$(grep -c INTEG inp-2); if [[ $m -eq 0 ]]; then nb=1; else nb=&
             &$(sed -n '/TROT/{n;p;}' inp-2); fi; geo=$(ls GEOMETRY*.xyz | hea&
             &d -1); [[ -z $geo ]] && exit 1; nat=$(wc -l $geo|cut -d' ' -f1);&
             &echo $nb $(($nat-2)) >> tmpin"
    msg = system(getcfg)
    if(msg /= 0) call stopgm(0, 'inp-2 not found in window ', dir(i))
    if(i == 1) then
        msg = system(getmol)
        if(msg /= 0) call stopgm(0, 'GEOMETRY*.xyz not found in window ', &
                                 dir(i))
    end if
    open(unit=11, file='tmpin')
    if (twham) then
        read(11, *, iostat=ioerr) ind, xi(i), ks(i)
    else
        read(11, *, iostat=ioerr) ind, xi(i)
    end if
    if(ioerr < 0) call stopgm("Incomplete inp-2 in dir ", dir(i))
    if(i == 1) then
        read(11, iostat=ioerr) nb, natom
        if(ioerr < 0) call stopgm("Incomplete inp-2 in dir ", dir(i))
    end if
    close(11, status='delete')
end subroutine

subroutine init_param()
#if defined(_OPENMP)
    use omp_lib
#endif
    use global
    use bluemoon
    use wham
    implicit none

    integer :: i, j, ind1, ind2
    real*8 :: x(nw)

    xi = xi * fac
    beta = 1 / kb / temp
    if(twham) then
        wbin = wbin * fac
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
        n = int((xmax - xmin) / wbin) ! total number of bins
        allocate(xbin(n))
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
    else
        invm = 1 / (invm * facm)
        do i = 1, 3
            do j = 1, 3
                invmij(i, j) = invm(i) * invm(j)
            end do
        end do
        kt = kb * temp
        allocate(dx(nw-1))
        do i = 1, nw-1
            dx(i) = abs((xi(i+1) - xi(i))) * 0.5
        end do
        n = nw
        allocate(xbin(n))
        xbin = xi
    end if
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
