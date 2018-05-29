subroutine stopgm(n, msgp1, msgp2)
    implicit none

    integer, intent(in) :: n
    character(len=*), intent(in) :: msgp1
    character(len=*), optional, intent(in) :: msgp2
    character(len=50) :: msg

    if(present(msgp2)) then
        msg = trim(msgp1 // adjustl(msgp2))
    else
        msg = msgp1
    end if 
    if(n > 0) call movefile(n, -1)
    write(*, *) 'WHAM stops because of ', msg
    stop
end subroutine
