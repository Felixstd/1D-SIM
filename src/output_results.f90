subroutine output_results(ts, expnb, solver, utp, zeta, eta)
  use size
  use resolution
  use global_var
  use shallow_water
  use MOMeqSW_output
  use rheology
  use option

  implicit none

  character filename*90

  integer :: i, k, Dt, Dx, adv
  integer, intent(in) :: ts, expnb, solver
  double precision, intent(in):: zeta(0:nx+1),eta(0:nx+1)
  double precision, intent(in)  :: utp(1:nx+1)
  double precision :: div(0:nx+1), sigma(0:nx+1), sig_norm(0:nx+1), zeta_norm(0:nx+1)
  double precision :: Erate(1:nx+1) ! KE rate loss/gain by rheology term

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
 elseif (adv_scheme .eq. 'semilag') then
    adv = 3
  endif

  Dt=int(Deltat) ! in s
  Dx=int(Deltax/1000d0) ! in km
  
  div(0) = 0d0
  div(nx+1) = 0d0
!  zeta(0) = 0d0
!  zeta(nx+1) = 0d0
  sigma(0) = 0d0
  sigma(nx+1) = 0d0
  sig_norm(0) = 0d0
  sig_norm(nx+1) = 0d0
  Erate = 0d0

  do i = 1, nx
     div(i) = (utp(i+1)-utp(i)) / Deltax ! calc divergence
     sigma(i) = (zeta(i)+eta(i))*div(i) - P_half(i)
     sig_norm(i) = (zeta(i)+eta(i))*div(i)*0.5d0/Pp_half(i) - 0.5d0*P_half(i)/Pp_half(i) ! norm by ice strength 
!     zeta_norm(i) = zeta(i) / (zmax_par*Pp_half(i))
  enddo
  
  do i = 2, nx
   Erate(i) = utp(i) * ( sigma(i) - sigma(i-1) ) / Deltax
  enddo
  
  write (filename, '("output/h_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') Dt, &
		    Dx,solver,IMEX, adv,BDF2,ts,expnb
  open (10, file = filename, status = 'unknown')
  
  write (filename, '("output/A_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (11, file = filename, status = 'unknown')

  write (filename, '("output/u_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (12, file = filename, status = 'unknown')

  write (filename, '("output/div_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (13, file = filename, status = 'unknown')

 write (filename, '("output/zeta_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
 open (14, file = filename, status = 'unknown')

  write (filename, '("output/eta_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
 open (15, file = filename, status = 'unknown')

!  write (filename, '("output/sigma_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (15, file = filename, status = 'unknown')

!  write (filename, '("output/zeta_norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (16, file = filename, status = 'unknown')

!  write (filename, '("output/sig_norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_ts",i4.4,".",i2.2)') Dt,Dx, &
!		    IMEX, adv,ts,expnb
!  open (17, file = filename, status = 'unknown')

  write (filename, '("output/Er_",i5.5,"s_",i5.5,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (18, file = filename, status = 'unknown')

  write (filename, '("output/Wdissip_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
 open (20, file = filename, status = 'unknown')

 if (mechenergy) then
  write (filename, '("output/Ppot_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (21, file = filename, status = 'unknown')
  
  write (filename, '("output/Pfric_R_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (22, file = filename, status = 'unknown')

  write (filename, '("output/Pfric_s_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (23, file = filename, status = 'unknown')

  write (filename, '("output/Plat_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (24, file = filename, status = 'unknown')

  write (filename, '("output/Ph_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (25, file = filename, status = 'unknown')
  
  write (filename, '("output/Pw_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i9.9,".",i2.2)') &
      Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (26, file = filename, status = 'unknown')
  endif


  write(10,*) ( h(i),       i = 0, nx+1 )
  write(11,*) ( A(i),       i = 0, nx+1 )
  write(13,*) ( div(i),     i = 0, nx+1 )
 write(14,*) ( zeta(i),    i = 0, nx+1 )
  write(15, *) (eta(i), i = 0, nx+1)
!  write(15,10) ( sigma(i),   i = 0, nx+1 )
!  write(16,10) ( zeta_norm(i),    i = 0, nx+1 )
!  write(17,10) ( sig_norm(i),   i = 0, nx+1 )
  write(12,*) ( utp(i),       i = 1, nx+1 )
  write(18,*) ( Erate(i),     i = 1, nx+1 )
  write(20,*) ( W_sigma(i),     i = 0, nx+1 )

  if (mechenergy) then
    write(21,*) ( P_pot(i),       i = 1, nx+1 )
    write(22,*) ( P_fric_R(i),     i = 1, nx+1 )
    write(23,*) ( P_fric_S(i),     i = 1, nx+1 )
    write(24,*) ( P_lat(i),     i = 1, nx+1 )
    write(25,*) ( P_h(i),     i = 1, nx+1 )
    write(26,*) ( P_w(i),     i = 1, nx+1 )
  endif

  do k = 10, 20
     close(k)
  enddo

  if (mechenergy) then 
    do k = 21, 26 
      close(k)
    enddo
  endif

  if (oceanSIM) then
  
  write (filename, '("output/etaw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (18, file = filename, status = 'unknown')

  write (filename, '("output/uw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (19, file = filename, status = 'unknown')
  
  write (filename, '("output/duwdt_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
		    Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (20, file = filename, status = 'unknown')
  
  write (filename, '("output/gdetawdx_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') & 
		    Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (21, file = filename, status = 'unknown')
  
  write (filename, '("output/tauaw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
		    Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (22, file = filename, status = 'unknown')
  
  write (filename, '("output/tauiw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
		    Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (23, file = filename, status = 'unknown')
  
  write (filename, '("output/buw_",i5.5,"s_",i3.3,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i6.6,".",i2.2)') &
		    Dt, Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (24, file = filename, status = 'unknown')
  
  
  write(18,10) ( etaw(i),      i = 0, nx+1 )
  write(19,10) ( uw(i),        i = 1, nx+1 )
  write(20,10) ( duwdt(i),     i = 1, nx+1 )
  write(21,10) ( gedetawdx(i), i = 1, nx+1 )
  write(22,10) ( tauaw(i),     i = 1, nx+1 )
  write(23,10) ( tauiw(i),     i = 1, nx+1 )
  write(24,10) ( buw(i),       i = 1, nx+1 )
  
  do k = 18, 24
     close(k)
  enddo
  
  endif
  
10 format (1x, 1000(f30.20, 1x))

  return
end subroutine output_results


subroutine output_sor(ts, k, solver, expnb, F)

  use size
  use resolution
  use global_var
  use shallow_water
  use MOMeqSW_output
  use rheology
  use option

    implicit none

    character filename*90


    integer :: i, Dt, Dx, adv
  integer, intent(in) :: ts, k, expnb, solver
  double precision, intent(in):: F(1:nx+1)
 
  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat) ! in min
  Dx=int(Deltax/1000d0) ! in km
  
  
  write (filename, '("output/b_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') Dt, &
		    Dx,solver, IMEX, adv,BDF2,ts,expnb
  open (11, file = filename, status = 'unknown')
  
  write(11,*) ( F(i), i = 1, nx+1 )

  close(11)
  
10 format (1x, 2000(f25.18, 1x))

  return

end subroutine

subroutine output_residual(ts, k, expnb, F)
  use size
  use resolution
  use option

  implicit none

  character filename*60

  integer :: i, Dt, Dx, adv
  integer, intent(in) :: ts, k, expnb
  double precision, intent(in):: F(1:nx+1)
 
  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km
  
  write (filename, '("output/res_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,"_ts",i4.4,"_k",i5.5,".",i2.2)') Dt, &
		    Dx,IMEX,adv,BDF2,ts,k,expnb
  open (11, file = filename, status = 'unknown')
  
  write(11,10) ( F(i), i = 1, nx+1 )

  close(11)
  
10 format (1x, 1000(f25.18, 1x))

  return
end subroutine output_residual

subroutine output_nb_ite(ts, k, fgmres_per_ts, expnb)
  use resolution
  use option
  implicit none

  character filename*60

  integer, intent(in) :: ts, k, fgmres_per_ts, expnb
  integer :: Dt, Dx, adv

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat) ! in s
  Dx=int(Deltax/1000d0) ! in km

  write (filename, '("output/Nbite_",i5.5,"s_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,"_ts",i4.4,".",i2.2)') Dt, &
	 Dx,IMEX,adv,BDF2,ts,expnb
  open (10, file = filename, access = 'append')
  
  write(10,10) ts,k-1,fgmres_per_ts

  close(10)

10 format (i5,1x,i4,1x,i5)

  return
end subroutine output_nb_ite

subroutine output_ini_L2norm(ts, L2norm, expnb)
  use resolution
  use option
  implicit none

  character filename*60

  integer :: Dt, Dx, adv
  integer, intent(in) :: ts, expnb
  double precision, intent(in) :: L2norm

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  write (filename, '("output/iniL2norm_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,"_ts",i4.4,".",i2.2)') Dt,&
	 Dx,IMEX,adv,BDF2,ts,expnb
  open (10, file = filename, access = 'append')
  
  write(10,10) ts,L2norm

  close(10)

10 format (i5,1x,f15.12)

  return
end subroutine output_ini_L2norm

subroutine output_u_and_du ( ts, k, utp, du )
  use size
  use resolution
!  use global_var
  use rheology
  use option
  implicit none

  character filename*60

  integer :: i, Dt, Dx, adv
  integer, intent(in) :: ts, k
  double precision, intent(in)  :: utp(1:nx+1), du(1:nx+1)

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  write (filename, '("output/uk1_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,"_ts",i4.4,"_k",i3.3,".dat")') Dt,&
		    Dx,IMEX, adv,BDF2,ts,k
  open (10, file = filename, status = 'unknown')
  
  write (filename, '("output/du_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1"_ts",i4.4,"_k",i3.3,".dat")') Dt,&
		    Dx,IMEX, adv,BDF2,ts,k
  open (11, file = filename, status = 'unknown')
  
  write(10,10) ( utp(i),       i = 1, nx+1 )
  write(11,10) ( du(i),       i = 1, nx+1 )

  close(10)
  close(11)

10 format (1x, 1000(f25.20, 1x))

  return
end subroutine output_u_and_du

subroutine output_diag_stress(ts, expnb, idiag)

! output diagnostic (ice-ocean and vice versa) stress at i=idiag

  use resolution
  use diag_stress
  use option
  implicit none

  character filename*70

  integer :: Dt, Dx, adv
  integer, intent(in) :: ts, expnb, idiag
  double precision :: ratio

  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
  endif

  Dt=int(Deltat/60d0) ! in min
  Dx=int(Deltax/1000d0) ! in km

  ratio = tauwidiag/tauiwdiag
  
  write (filename, '("output/diag_stress_",i3.3,"min_",i3.3,"km_IMEX",i1.1,"_adv",i1.1,"_BDF2_",i1.1,".",i2.2)') Dt,&
	 Dx,IMEX,adv,BDF2,expnb
  open (10, file = filename, access = 'append')
  
  write(10,10) ts, idiag, tauwidiag, tauiwdiag, ratio, tauaidiag

  close(10)

10 format (i5,1x, i5,1x,f12.8,1x,f12.8,1x,f12.8,1x,f12.8)

  return
end subroutine output_diag_stress


subroutine output_times(ts, expnb, solver, nstep, A_time, h_time, u_time)

  use size
  use resolution
  use global_var
  use shallow_water
  use MOMeqSW_output
  use rheology
  use option

  implicit none

  character filename*90

  integer :: i, k, Dt, Dx, adv
  integer, intent(in) :: ts, expnb, solver, nstep
  double precision, intent(in) :: A_time(:), h_time(:), u_time(:)


  if (adv_scheme .eq. 'upwind') then
    adv = 1
  elseif (adv_scheme .eq. 'upwindRK2') then
    adv = 2
 elseif (adv_scheme .eq. 'semilag') then
    adv = 3
  endif

  Dt=int(Deltat) ! in s
  Dx=int(Deltax/1000d0) ! in km

    
  write (filename, '("output/h_time_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
        Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
  open (10, file = filename, status = 'unknown')

  write (filename, '("output/A_time_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
		    Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
  open (11, file = filename, status = 'unknown')

    write (filename, '("output/u_time_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
		    Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
  open (12, file = filename, status = 'unknown')

  write(10,*) ( A_time(i),       i = 1, nstep+1 )
  write(11,*) ( h_time(i),       i = 1, nstep+1 )
  write(12,*) ( u_time(i),       i = 1, nstep+1 )

  close(10)
  close(11)
  close(12)

end subroutine output_times


subroutine output_mech_energy(P_pot_Tavg, P_fric_R_Tavg, P_fric_S_Tavg, &
                              P_fric_R_visc_Tavg, P_fric_S_visc_Tavg, P_fric_R_plas_Tavg,&
                             P_fric_S_plas_Tavg, P_lat_Tavg, P_h_Tavg, P_w_Tavg, &
                             ts, expnb, solver, nstep)
    
    use size
    use resolution
    use global_var
    use rheology
    use option

    implicit none 
    character filename*90

    integer :: i, k, Dt, Dx, adv, tstep
    integer, intent(in) :: ts, expnb, solver, nstep
    double precision, intent(in) :: P_pot_Tavg(nstep), P_fric_R_Tavg(nstep), P_fric_S_Tavg(nstep)
    double precision, intent(in) :: P_fric_R_visc_Tavg(nstep), P_fric_S_visc_Tavg(nstep)
    double precision, intent(in) :: P_fric_R_plas_Tavg(nstep), P_fric_S_plas_Tavg(nstep)
    double precision, intent(in) :: P_lat_Tavg(nstep), P_h_Tavg(nstep), P_w_Tavg(nstep)


    if (adv_scheme .eq. 'upwind') then
        adv = 1
    elseif (adv_scheme .eq. 'upwindRK2') then
        adv = 2
    elseif (adv_scheme .eq. 'semilag') then
        adv = 3
    endif

    Dt=int(Deltat) ! in s
    Dx=int(Deltax/1000d0) ! in km
    

        
    write (filename, '("output/Pp_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
            Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (10, file = filename, status = 'unknown')

    write (filename, '("output/PfR_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (11, file = filename, status = 'unknown')

    write (filename, '("output/PfS_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (12, file = filename, status = 'unknown')

    write (filename, '("output/Pl_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
            Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (13, file = filename, status = 'unknown')

    write (filename, '("output/Ph_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (14, file = filename, status = 'unknown')

    write (filename, '("output/Pw_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (15, file = filename, status = 'unknown')

    write (filename, '("output/PfRp_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (16, file = filename, status = 'unknown')

    write (filename, '("output/PfSp_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (17, file = filename, status = 'unknown')

    write (filename, '("output/PfRv_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (18, file = filename, status = 'unknown')

    write (filename, '("output/PfSv_t_",i5.5,"s_",i6.6,"km_solv",i1.1,"_IMEX",i1.1,"_adv",i1.1,"_BDF2",i1.1,"_ts",i8.8,".",i2.2)') &
                Dt, Dx,solver,IMEX, adv,BDF2,ts,expnb
    open (19, file = filename, status = 'unknown')


    write(10,*) ( P_pot_Tavg(i),      i = 1, nstep )
    write(11,*) ( P_fric_R_Tavg(i),   i = 1, nstep)
    write(12,*) ( P_fric_S_Tavg(i),   i = 1, nstep )
    write(13,*) ( P_lat_Tavg(i),      i = 1, nstep )
    write(14,*) ( P_h_Tavg(i),        i = 1, nstep)
    write(15,*) ( P_w_Tavg(i),        i = 1, nstep )
    write(16,*) ( P_fric_R_plas_Tavg(i),        i = 1, nstep )
    write(17,*) ( P_fric_S_plas_Tavg(i),        i = 1, nstep )
    write(18,*) ( P_fric_R_visc_Tavg(i),        i = 1, nstep )
    write(19,*) ( P_fric_S_visc_Tavg(i),        i = 1, nstep )


    do i = 10, 19
        close(i)
    enddo


end subroutine