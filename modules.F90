module global
    implicit none
    integer :: n, nw, nb, natom, ncut, nskip, nread, ind(3)
    real*8, parameter :: kb = 0.316679d-5 ! boltzmann constant in au/kelvin
    real*8, parameter :: fac = 1.889725989 ! from angstrom to au
    real*8, parameter :: fac1 = 6.27509d2 ! from au to kcal/mol
    real*8, allocatable :: xi(:), xbin(:)
    real*8 :: temp, beta, kt
    character(len=200) :: rootdir, date*30
    character(len=20), allocatable :: dir(:)
    logical :: symm, tbl, twham, tdiff, tproj
    contains
    subroutine stopgm(msgp1, msgp2)
        implicit none
    
        character(len=*), intent(in) :: msgp1
        character(len=*), optional, intent(in) :: msgp2
        character(len=50) :: msg
    
        if(present(msgp2)) then
            write(msg, '(3a)') msgp1, msgp2, '.'
        else
            write(msg, '(2a)') msgp1, '.'
        end if
        ! for some unknown reasons, ifort does not stop + arg
#if defined(__INTEL_COMPILER)
        write(*, '(a)') msg
        stop
#else
        stop msg
#endif
    end subroutine
end module global

module bluemoon
    implicit none
    real*8, allocatable :: dx(:)
    real*8 :: invm(3), invmij(3, 3)
    real*8, parameter :: h = 0.001, h2 = 0.5 / h , h42 = 0.25 / h ** 2
end module

module wham
    implicit none
    real*8, parameter :: facm = 1822.888 ! convert the unit of mass
    real*8, allocatable :: ni(:), ks(:)
    real*8 :: wbin, xmax, xmin
    real*8 :: tol
end module
