program test_read
    implicit none
    integer i(4), clock

    open(unit=11, file='test')
    call test_unit()
end program

subroutine test_unit()
    implicit none

    integer i(4), clock
    read(11, *) i, clock
    write(*, *) i, clock
end subroutine
