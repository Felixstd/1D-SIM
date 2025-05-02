subroutine ini_get(utp)

    use size
    use global_var
    use muphi
    use option

    implicit none 

    double precision, intent(inout)  :: utp(0:nx+1)
    integer :: i

    utp(1)    = 0d0 ! close bc
    utp(nx)   = 0d0
    utp(nx+1) = 0d0 ! close bc
    h(0)    = 0d0
    h(nx+1) = 0d0
    A(0)    = 0d0
    A(nx+1) = 0d0

    do i = 1, nx
        !10000
        A(i) = min(max((real(40000d0)-(real(i)-real(nx)/2d0)**2)**(1/20d0), 0d0), 1d0)

        if (A(i) .gt. 0d0) then  
            h(i) = min(max((real(40000d0)-(real(i)-real(nx)/2d0)**2)**(1/20d0), 0d0), 1d0)
            ! print*, real(i)
            utp(i) = real(i)/1000d0
            print*, utp(i)
        endif

        ! if (h(i) > 0d0) then
        !     utp(i) =  0.01d0!min(max((real(10000d0)-(real(i)-real(nx)/2d0)**2)**(1/20d0), 0d0), 1d0)/100d0
        ! endif
        ! utp(i) = 10d0
    enddo   

    utp(1) = 0d0
    utp(nx) = 0d0
    utp(nx+1) = 0d0



end subroutine ini_get