subroutine stopgm(n, msg)
    implicit none

    integer, intent(in) :: n
    character, intent(in) :: msg

    if(n > 0) call movefile(n, -1)
    write(*, *) 'WHAM stops because of ', msg
    stop
end subroutine
