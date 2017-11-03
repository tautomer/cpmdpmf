subroutine read_conf()
    use global
    implicit none

    open(unit=20,file="config", action="read")
    read(20,*) wbin ! unit in angstrom
    read(20,*) ncut
    read(20,*) ks
    read(20,*) temp
    read(20,*) nw
    read(20,*) tol
    read(20, *) dir
    close(20)

end subroutine

subroutine read_meta()
    use global
    implicit none

    integer i
    character(len=12) :: path, rmtmp
    character(len=70) :: getcfg
    character(len=200) :: getmol

    getcfg = "grep DIF inp-2 | cut -d' ' -f4-7 > tmpin; sed -n '/MAXS/{n;p;}'&
             & inp-2 >> tmpin"
    getmol = "m=$(grep -c TROT inp-2); if [[ $m -eq 0 ]]; then nb=1; else nb=$(&
             &sed -n '/TROT/{n;p;}' inp-2); fi; nat=$(wc -l $(ls GEOMETRY*.xyz &
             &| head -1) | cut -d' ' -f1); echo $nb $nat >> tmpin"
    rmtmp = "rm -f tmpin"
    call getcwd(rootdir)
    do i = 1, nw
        path = trim("./" // adjustl(dir(i)))
        call chdir(path)
        call system(getcfg)
        if(i.eq.1) call system(getmol)
        open(unit=11, file='tmpin')
        read(11, *) ind, xi(i)
        read(11, *) nsteps(i)
        if(i.eq.1) nb, natom
        close(11)
        call chdir(rootdir)
    end do
    n = dint((xmax - xmin) / wbin) ! total number of bins
    wbin = wbin * fac
    xmax = (3 * xi(nw) - xi(nw-1)) * fac * 0.5
    xmin = (3 * xi(1) - xi(2)) * fac * 0.5
    beta = 1.d0 / kb / temp
    !!! need a folder to store data of each window
end subroutine
