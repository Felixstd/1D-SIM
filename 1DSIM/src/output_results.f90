subroutine output_results(ts, expnb, utp, zeta, eta)

    use size
    use resolution
    use global_var
    use rheology


    implicit none 

    character filename*82


    integer :: i, k, Dt, Dx
    integer, intent(in) :: ts, expnb
    double precision, intent(in):: zeta(0:nx+1),eta(0:nx+1)
    double precision, intent(in)  :: utp(0:nx+1)

    Dt=int(Deltat) ! in s
    Dx=int(Deltax/1000d0) ! in km

    write (filename, '("output/h_",i5.5,"km",i5.5,"_ts",i6.6,".",i2.2)') Dt, &
            Dx,ts,expnb
    open (10, file = filename, status = 'unknown')

    write (filename, '("output/A_",i5.5,"km",i5.5,"_ts",i6.6,".",i2.2)') Dt, &
            Dx,ts,expnb
    open (11, file = filename, status = 'unknown')

    write (filename, '("output/u_",i5.5,"km",i5.5,"_ts",i6.6,".",i2.2)') Dt, &
            Dx,ts,expnb
    open (12, file = filename, status = 'unknown')

    write(10,*) ( h(i),       i = 0, nx+1 )
    write(11,*) ( A(i),       i = 0, nx+1 )
    write(12,*) ( utp(i),     i = 0, nx+1 )

    do k = 10, 12
        close(k)
    enddo

10 format (1x, 1000(f25.20, 1x))
end subroutine output_results