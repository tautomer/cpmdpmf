module global
    integer n, nw, nb, natom
    integer ncut
    integer, allocatable :: nsteps(:)
    real*8, parameter :: kb = 0.316679d-5 ! boltzmann constant in au/kelvin
    real*8, parameter :: fac = 1.889725989 ! from angstrom to au
    real*8, parameter :: fac1 = 627.509d0 ! from au to kcal/mol
    real*8, allocatable :: xi(:), ni(:), xbin(:)
    real*8 ks, wbin, temp, beta, tol, ind(3)
    real*8 xmax, xmin
    character(len=100) :: rootdir
    character(len=10), allocatable :: dir(:)
end module global
