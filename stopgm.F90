subroutine stopgm(msgp1, msgp2)
    implicit none

    character(len=*), intent(in) :: msgp1
    character(len=*), optional, intent(in) :: msgp2
    character(len=50) :: msg

    if(present(msgp2)) then
        write(msg, '(3a)') msgp1, msgp2, '.'
    else
        msg = msgp1
    end if
    ! for some unknown reasons, ifort does not stop + arg
#if defined(__INTEL_COMPILER)
    write(*, '(a)') msg
    stop
#else
    stop msg
#endif
end subroutine
