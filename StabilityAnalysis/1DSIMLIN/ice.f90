program ice

    use size 
    use rheology
    use muphi
    use properties
    use resolution
    use global_var
    use numerical
    use option

    implicit none 


    logical :: restart
    integer :: i, ii, ts, tsini, nstep, tsfin, k, s, Nmax_OL, solver, idiag
    integer :: out_step(11), expnb, expres, ts_res
    integer :: time_out, n_output, output_time
    double precision :: e
    double precision :: u(0:nx+1), un1(0:nx+1), un2(0:nx+1)
    double precision :: zeta(0:nx+1), eta(0:nx+1), mn1(0:nx+1), m(0:nx+1)
    double precision :: meanvalue, time1, time2, timecrap
    double precision :: nbhr

    !------ Setting Constants -------!
    expnb = 1
    restart = .false.
    solver = 'explicit'

    T_tot = 1*24*60*60
    Deltat = 1
    nstep = T_tot/Deltat

    write (filename,'("output/time_run.",i2.2)') expnb
    open (10, file = filename, status = 'unknown')
    time_out = 0
    output_time = 0
    write(10, *) (output_time)
    out_step(1) = 0
    do n_output = 2, 11
        time_out = time_out + T_tot / 10
        output_time = output_time + nstep/10
        write(10, *) (output_time)
        out_step(n_output)=output_time   
    enddo
    close(10)


end program ice 