subroutine viscouscoefficient(utp, zeta, eta)
	use size
	use resolution
	use rheology
	use option
	use global_var
	use muphi
    use properties
	
	implicit none

	integer :: i

	double precision, intent(in) :: utp(1:nx+1)
	double precision, intent(out):: zeta(0:nx+1), eta(0:nx+1)
	double precision :: dudx, deno, denonum, denomin
	double precision :: shearmax, Inertial_num, muI, zeta_i, eta_i 
! 
	! denomin=2d-09
  	denomin = 2d-05

	muI = 0d0
	Inertial_num = 0d0


	if (rheo .eq. 2) then

		do i = 1, nx
			zeta_i =0 
			eta_i = 0
			dudx = ( utp(i+1) - utp(i) ) / Deltax
			shearmax = sqrt( (dudx)**2d0+small2) 

			Inertial_num = min(d_average * shearmax * SQRT(rho * h(i)/(Pp_half(i)+Tp_half(i))), 1d0)
			muI = max(mu_0 + ( mu_infty - mu_0 ) / ( I_0/Inertial_num + 1), mu_0)

			! muI = mu_0 + mu_infty*Inertial_num
			mu_I(i) = muI

			if (regularization .eq. 'tanh') then

				zeta(i) = 2d08*Pp_half(i) * tanh((mu_b) / (shearmax * 2d08))
				eta(i)  = eta_max * tanh((mu_I(i) / 2) * Pp_half(i) / (shearmax * (eta_max)))

			elseif (regularization .eq. 'capping') then
				
				zeta(i) = min(zeta_max,  (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax))

				eta(i) = min(eta_max,  (muI) *(Pp_half(i)+Tp_half(i)) / (2d0*shearmax))

			elseif (regularization .eq. 'ringeisen') then 

				zeta_i = (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax)
				eta_i = (muI) *(Pp_half(i)+Tp_half(i)) / (2d0*shearmax)
				
				zeta(i) = min(zeta_i, &
						zeta_max*min(1d0,zeta_i/eta_i))
				eta(i) = min(eta_i, &
						eta_max*min(1d0,eta_i/zeta_i))

			elseif (regularization .eq. 'ringeisen2') then 
				zeta_i = (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax)
				eta_i = (muI) *(Pp_half(i)+Tp_half(i)) / (2d0*shearmax)
				
				zeta(i) = min(zeta_max,  (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax))
				eta(i)  = min(eta_i, zeta_max*eta_i/zeta_i, eta_max)
			
			elseif (regularization .eq. 'div') then 

				zeta_i = (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax)
				eta_i = (muI) *(Pp_half(i)+Tp_half(i)) / (2d0*shearmax)

				zeta(i) = min(zeta_i, zeta_max* min(1d0, eta_i/eta_max, zeta_i/zeta_max))
				eta(i) = min(eta_i, eta_max* min(1d0, eta_i/eta_max, zeta_i/zeta_max))
				! if (shearmax < 1d-10) then 
				! 	zeta(i) = zeta_max
				! 	eta(i)  = eta_max
				! else
				! 	zeta(i) = (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax)
				! 	eta(i) = (muI) *(Pp_half(i)+Tp_half(i)) / (2d0*shearmax)
				! endif
			
			elseif (regularization .eq. 'hard') then 
				zeta(i) = min(zeta_max,  (mu_b) *(Pp_half(i)+Tp_half(i)) / (shearmax))

				if (zeta(i) .eq. zeta_max) then
					eta(i) = eta_max
				else 
					eta(i) = (muI) *(Pp_half(i)+Tp_half(i)) / (2d0*shearmax)

				endif
			
			endif

			if (rep_closure) then 
				P_half(i) = 2*abs(dudx)*zeta(i)
			else


				P_half(i) = Pp_half(i)-Tp_half(i)
			endif


		enddo

	else 

		do i = 1, nx ! for tracer points

			dudx = ( utp(i+1) - utp(i) ) / Deltax

			if ( regularization .eq. 'tanh' ) then

				deno = alpha*sqrt( (dudx)**2d0 + small2 ) ! small2 is there to avoid div by zero
				!     deno = alpha*(abs(dudx))
				!     deno = max( deno, 1d-30 )
				! if (Pstart_change) then 

				zeta(i) = ((Pp_half(i)+Tp_half(i))/denomin)*tanh(denomin*(1d0/deno))

			elseif ( regularization .eq. 'Kreyscher' ) then

				deno = alpha*sqrt( (dudx)**2d0)
				zeta(i) = (Pp_half(i)+Tp_half(i)) / ( deno + denomin )

			elseif ( regularization .eq. 'capping' ) then
      

        		zeta(i) = min(Pp_half(i)/(alpha*sqrt( (dudx)**2d0+small2)), zeta_max)

				! if (Pstart_change) then 
				! 	zeta(i) = min(Pp_half(i)/(alpha*sqrt( (dudx)**2d0+small2)), zeta_max)
				! endif
				! print*, zeta(i), zeta(i) * e_2, Pp_half(i), sqrt( (dudx)**2d0+small2), alpha
        		! print*, zeta(i)
			
			elseif (regularization .eq. 'viscous') then 
				zeta(i) = zeta_max


      		else
        		print *, 'Wrong regularization'
				stop
        
      		endif

      		eta(i)  = zeta(i) * e_2

			if (rep_closure) then  ! replacement closure (Kreysher et al. 2000)
		
				if ( regularization .eq. 'tanh' ) then
					P_half(i) = (Pp_half(i)-Tp_half(i)) * ( deno / denomin ) * tanh(denomin*(1d0/deno))
			
				elseif ( regularization .eq. 'Kreyscher' ) then
					P_half(i) = (Pp_half(i)-Tp_half(i)) * deno / ( deno + denomin )
			
				elseif ( regularization .eq. 'capping' ) then
					denonum=alpha*sqrt( (dudx)**2d0)
					P_half(i) = (Pp_half(i)-Tp_half(i)) * denonum / deno
			
				endif
        
			else
		
				P_half(i) = Pp_half(i)-Tp_half(i) ! P_half includes tensile strength
    
      		endif

    	enddo

    	if (linear_viscous) then
      		zeta = 1d08
      		eta  = e_2 * zeta
      		P_half = 0d0
    	endif
	
	endif 

	return
end subroutine viscousCoefficient



