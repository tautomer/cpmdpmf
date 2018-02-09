module global
    integer n, nw, nb, natom, ncut, nskip, ind(3)
    integer, allocatable :: nsteps(:)
    real*8, parameter :: kb = 0.316679d-5 ! boltzmann constant in au/kelvin
    real*8, parameter :: fac = 1.889725989 ! from angstrom to au
    real*8, parameter :: fac1 = 6.27509d2 ! from au to kcal/mol
    real*8, allocatable :: xi(:), ni(:), xbin(:), ks(:)
    real*8 wbin, temp, beta, tol
    real*8 xmax, xmin
    character(len=200) :: rootdir
    character(len=10), allocatable :: dir(:)
end module global
