subroutine laplacian(x, lap_x)

    use size
    use resolution
    use properties
    use numerical
    use global_var

    implicit none

    integer :: i
    double precision, intent(in) :: x(1:nx+1)
    double precision, intent(out):: lap_x(1:nx+1)

    lap_x(1)    = 0d0 ! close bc
    lap_x(nx+1) = 0d0 ! close bc 

    do i = 2, nx

        lap_x(i) = (x(i+1) - 2d0*x(i) + x(i-1))/(Deltax2)

        if (lap_x(i) < 1d-15) lap_x(i) = 0d0
        
    enddo

end subroutine