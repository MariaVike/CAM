module gw_chem

!
! This module contains code to compute the wave‐driven constituent transport 
! due to non-breaking gravity waves (Gardner et al, ESS, 2019 ; Gardner et al, JGR, 2018).
! Here we compute the effective wave diffusivity (K_wave) as a function of the 
! eddy diffusivity (Kzz) and of thevariances of the temperature and lapse rate fluctuations
!

use gw_utils,     only: r8
use coords_1d,    only: Coords1D
use spmd_utils,   only: masterproc
use cam_logfile,  only: iulog

implicit none
private
save

! Public interface.
public :: effective_gw_diffusivity

contains

!==========================================================================

subroutine effective_gw_diffusivity (ncol, band, lambda_h, p, dt,    &
           t, rhoi, nm, ni, c_speed, tau, egwdffi, ubi, k_wave, xi, k_e, k_eff,  &
           zm, zi, var_t, dtdz, lat, kwvrdg)
!-----------------------------------------------------------------------
! Compute K_wave (wave effective diffusivity) etc...
! ....
!-----------------------------------------------------------------------
use gw_common, only: GWBand, pver, pi
use physconst, only: cpair, cpairv, gravit
use gw_utils, only: midpoint_interp
!------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! band of phase speeds c within whic waves are emitted
  type(GWBand), intent(in) :: band
  ! horizonatl wavelength
  real(r8), intent(in) :: lambda_h
  ! horiz wavenumber [anisotropic orography].
  real(r8), intent(in), optional :: kwvrdg(ncol)
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Time step.
  real(r8), intent(in) :: dt
  ! Midpoint temperature.
  real(r8), intent(in) :: t(ncol,pver)
  ! Interface density (Kg m-3)
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisala frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at interfaces.
  real(r8), intent(in) :: ubi(ncol,pver+1)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c_speed(ncol,-band%ngwv:band%ngwv)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(in) :: egwdffi(ncol,pver+1)
  ! Midpoint and Interface altitudes above ground (m).
  real(r8), intent(in) :: zm(ncol,pver), zi(ncol,pver+1)    
  ! Latitude in radians.
  real(r8), intent(in) :: lat(:)


  real(r8), intent(out) :: k_wave(ncol,pver) !total over entire wave spectrum for each GW source (i.e. Beres and C&M)
  real(r8), intent(out) :: xi(ncol,pver)
  real(r8), intent(out) :: k_e(ncol,pver)
  real(r8), intent(out) :: k_eff(ncol,pver)

  ! Variance of T'
  real(r8), intent(out)  :: var_t(ncol,pver)
  ! Dt/Dz environment 
  real(r8), intent(out)  :: dtdz(ncol, pver)

  
  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, and column loop indices.
  integer  :: k, l, i
  real(r8) :: icount, times
  ! Gas consant for dry air (m2 K-1 s-2)
  real(r8), parameter :: R_air = 287._r8
  ! Adiabatic lapse rate and R/cp ratio
  real(r8) :: gamma_ad, cp_r 
  ! Interface temperature.
  real(r8) :: ti(ncol,pver+1)
  ! The absolute momenum flux computed from tau
  real(r8) :: mom_flux(ncol,-band%ngwv:band%ngwv,pver+1)
  ! GW intrinsic frequency 
  real(r8):: gw_frq(ncol,-band%ngwv:band%ngwv, pver+1)
  ! GW intrinsic phase speed
  real(r8):: c_i(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Vertical wavenumber m
  real(r8):: m(ncol,-band%ngwv:band%ngwv,pver+1) 
  !convert degrees to radians
  real(r8), parameter :: degree_radian = pi/180._r8
  !convert hours to seconds
  real(r8), parameter :: hr_sec = 3600._r8
  !Coriolis frequency 
  real(r8) :: coriolis_f(size(lat))

  ! Temporary values for wave filtering
  real(r8):: var_t_initial(ncol,pver+1)

  ! Temporary values used for calculations
  real(r8), dimension(ncol, pver)  :: f_n_gammad, f_ln_nf,lapse_rate_sq, &
				      k_eff_sqr, one_min_xi

 !compute adiabatic lapse rate(K/m) and R/Cp ratio
 gamma_ad=gravit/cpair
 cp_r= cpair/R_air
 !compute temperature at interface (call function)
 ti(:,2:pver)=midpoint_interp(t)

 !Compute the coriolis frequency (rad/s) and set it to 2pi/24h for equatorial regions
  where (abs(lat) <= 30._r8*degree_radian)
   coriolis_f=(2._r8*pi)/(24._r8*hr_sec)
  elsewhere
   coriolis_f=abs( (2._r8*pi*sin(lat))/(12._r8*hr_sec) )
  end where
 
 
do i=1,ncol
 do k = 2, pver
  do l = -band%ngwv, band%ngwv
        !compute momentum flux MF (m2/s2) from the wave stress tau (Pa)
         mom_flux(i,l,k)=tau(i,l,k)/rhoi(i,k)     

	 !compute gw intrinsic frequency and vertical wavenumber
         c_i(i,l,k)=c_speed(i,l)-ubi(i,k) 

         IF (c_i(i,l,k) .ne. 0) then
	   gw_frq(i,l,k)=abs(c_i(i,l,k))/lambda_h

           if (present(kwvrdg)) then
            m(i,l,k)= (ni(i,k)*kwvrdg(i))/gw_frq(i,l,k)
	   else
            m(i,l,k)= (ni(i,k)*band%kwv)/gw_frq(i,l,k)
	   endif
         ELSE
          !at critical levels where c_i=0 (i.e. c=ubi) tau=0
	  !set everything else to zero too otherwise lateron in the
	  !computation of m we get division by zeros and generation of NaNs and Inf
          gw_frq(i,l,k)=0._r8  
	  m(i,l,k)=0._r8
         ENDIF
      
  enddo
 enddo
enddo

!Compute the Total temperature perturbation variance and the total effective diffusivity for the initial spectrum 
 IF (present(kwvrdg)) then          
  call compute_VarT_Keff (ncol, band, lambda_h, t, ti, rhoi, ni, egwdffi, &
           		  zm, cp_r, gamma_ad, coriolis_f,  mom_flux,      &
           	 	  gw_frq, m, var_t, dtdz, k_eff, k_eff_sqr,       &
           		  lapse_rate_sq, f_n_gammad, f_ln_nf, kwvrdg=kwvrdg)
 ELSE
   call compute_VarT_Keff (ncol, band, lambda_h, t, ti, rhoi, ni, egwdffi, &
           		  zm, cp_r, gamma_ad, coriolis_f,  mom_flux,      &
           	 	  gw_frq, m, var_t, dtdz, k_eff, k_eff_sqr,       &
          		  lapse_rate_sq, f_n_gammad, f_ln_nf)
 ENDIF 

!times=0.
wave_filter: DO
!times=times+1
 var_t_initial(:,:)=var_t(:,:)
 do i=1,ncol
  do k = 2, pver
   do l = -band%ngwv, band%ngwv
     !filter out from spectrum those waves that are damped by wave-induced diffusion and retain only
     !those satisfying the following criterion 
     if ( m(i,l,k)*k_eff(i,k) .lt. gw_frq(i,l,k)/m(i,l,k)) then 
        m(i,l,k)=m(i,l,k)
        gw_frq(i,l,k)=gw_frq(i,l,k)
     else
        m(i,l,k)=0._r8
        gw_frq(i,l,k)=0._8
     endif
    enddo
  enddo
 enddo


 !Recompute var_t and k_eff for the new filtered spectrum
IF (present(kwvrdg)) then          
  call compute_VarT_Keff (ncol, band, lambda_h, t, ti, rhoi, ni, egwdffi, &
           		  zm, cp_r, gamma_ad, coriolis_f,  mom_flux,      &
           	 	  gw_frq, m, var_t, dtdz, k_eff, k_eff_sqr,       &
           		  lapse_rate_sq, f_n_gammad, f_ln_nf, kwvrdg=kwvrdg)
 ELSE
   call compute_VarT_Keff (ncol, band, lambda_h, t, ti, rhoi, ni, egwdffi, &
           		  zm, cp_r, gamma_ad, coriolis_f,  mom_flux,      &
           	 	  gw_frq, m, var_t, dtdz, k_eff, k_eff_sqr,       &
          		  lapse_rate_sq, f_n_gammad, f_ln_nf)
 ENDIF

  !evaluate if the new var_t is much different from the previous var_t, IF change is larger than 5% of previous 
  !value anywhere in the domain (icount not zero), re-filter wave spectrum (loop continues). If change is less than 5%
  !everywhere (icount=0) var_t and k_eff stabilized thus exit loop.
  icount=0.
  do i=1,ncol
   do k = 2, pver
     if (  (abs(var_t(i,k)-var_t_initial(i,k)) .gt. 0.05*var_t_initial(i,k)) ) then  
        icount=icount+1												 
     endif
   !IF (masterproc) then
   !    if (k .eq. 14) then
   !     if (var_t_initial(i,k) .ne. 0._r8) then
   !    write (iulog,*) 'var_t, var_t_initial, k_eff', i, k, var_t(i,k), var_t_initial(i,k), k_eff(i,k)
   !    write (iulog,*) 'var_t - var_t_initial', var_t(i,k)- var_t_initial(i,k)
   !    write (iulog,*) 'times, icount', times, icount
   !    endif
   !    endif
   !ENDIF
  enddo
 enddo

  if (icount .eq. 0.) then
       !IF (masterproc) then
        !write (iulog,*) 'EXITING WAVE FILTER HERE icount=', icount, times
       !ENDIF
    exit wave_filter
  endif

ENDDO wave_filter

 !compute instability parameter and energy term using the final Var(T) and k_eff values after filtering
 do k=2,pver
     where (k_eff(:,k) .ne. 0.) 
      k_e(:,k)=(var_t(:,k)/lapse_rate_sq(:,k))*( (4._r8*f_n_gammad(:,k))/ti(:,k) )*k_eff_sqr(:,k)*(coriolis_f)**0.5
      xi(:,k)=(var_t(:,k)/lapse_rate_sq(:,k))*( f_ln_nf(:,k)/(2*k_eff(:,k)) ) 
      xi(:,k)=min(0.99, xi(:,k))
     elsewhere
      k_e(:,k)=0.
      xi(:,k)=0.
     end where
  enddo


 !Finally compute k_wave
   do k = 2, pver
      one_min_xi(:,k) =1._r8-xi(:,k)    
      k_wave(:,k)=(1._r8/one_min_xi(:,k))*( xi(:,k)*egwdffi(:,k)+ & 
                   (cp_r-1._r8)*k_e(:,k) )
   enddo

   !set variables at model top (k=1) to be zero
   k_wave(:,1)=0._r8
   xi(:,1)=0._r8
   k_e(:,1)=0._r8
   k_eff(:,1)=0._r8
   !k_wave(:,pver)=0._r8
   !xi(:,pver)=0._r8
   !k_e(:,pver)=0._r8
   !k_eff(:,pver)=0._r8


end subroutine effective_gw_diffusivity

!==========================================================================
subroutine compute_VarT_Keff (ncol, band, lambda_h, t, ti, rhoi, ni, egwdffi, &
           		      zm, cp_r, gamma_ad, coriolis_f,  mom_flux,      &
           	 	      gw_frq, m, var_t, dtdz, k_eff, k_eff_sqr,       &
           		      lapse_rate_sq, f_n_gammad, f_ln_nf, kwvrdg)

!-----------------------------------------------------------------------
! Compute Var(T') and K_eff ...
! ....
!-----------------------------------------------------------------------
use gw_common, only: GWBand, pver, pi
use physconst, only: gravit
!------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! band of phase speeds c within wich waves are emitted
  type(GWBand), intent(in) :: band
  ! horizonatl wavelength
  real(r8), intent(in) :: lambda_h
  ! Midpoint temperature.
  real(r8), intent(in) :: t(ncol,pver)
  ! Interface temperature.
  real(r8), intent(in) :: ti(ncol,pver+1)
  ! Interface density (Kg m-3)
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Interface Brunt-Vaisala frequency.
  real(r8), intent(in) ::  ni(ncol,pver+1)
  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(in) :: egwdffi(ncol,pver+1)
  ! Midpoint and Interface altitudes above ground (m).
  real(r8), intent(in) :: zm(ncol,pver)   
  !Coriolis frequency 
  real(r8), intent(in)  :: coriolis_f(:)
  ! Adiabatic lapse rate and R/cp ratio
  real(r8), intent(in) :: gamma_ad, cp_r 
  ! The absolute momenum flux computed from tau
  real(r8), intent(in) :: mom_flux(ncol,-band%ngwv:band%ngwv,pver+1)
  ! GW intrinsic frequency 
  real(r8), intent(in) :: gw_frq(ncol,-band%ngwv:band%ngwv, pver+1)
  ! Vertical wavenumber m
  real(r8), intent(in) :: m(ncol,-band%ngwv:band%ngwv,pver+1) 
  ! horiz wavenumber [anisotropic orography].
  real(r8), intent(in), optional :: kwvrdg(ncol)
  

  real(r8), intent(out) :: k_eff(ncol,pver)
  real(r8), intent(out) :: k_eff_sqr(ncol,pver)
  real(r8), intent(out) :: var_t(ncol,pver)
  real(r8), intent(out) :: dtdz(ncol, pver)  
  ! Temporary values used for calculations also to be used for the computaton of k_e and xi
  real(r8), intent(out) :: lapse_rate_sq(ncol,pver)
  real(r8), intent(out) :: f_n_gammad(ncol,pver)
  real(r8), intent(out) :: f_ln_nf(ncol,pver)

  
  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  integer :: k, l, i
  ! The vertical wavelength and the lmbd_h/lmbd_z ratio
  real(r8) :: lambda_z(ncol,-band%ngwv:band%ngwv,pver+1) 
  real(r8) :: lambda_ratio(ncol,-band%ngwv:band%ngwv,pver+1) 
  ! horiz wavelength  [anisotropic orography]
  real(r8) :: lambda_h_rdg(ncol) 
  ! GW temperature perturbation
  real(r8):: gw_t(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Temporary values used for calculations locally
  real(r8), dimension(ncol,pver+1) :: g_NT_sq
  real(r8), dimension(ncol, pver)  :: b, b_sq, c
				      
 
!compute variance of tot temperature perturbation Var(T') = SUM(T'^2) across whole spectrum 
do i=1,ncol
 do k = 2, pver
  var_t(i,k)=0._r8
  do l = -band%ngwv, band%ngwv

      ! FOR INITIAL SPECTRUM: exclude critical levels where c_i=0 and thus gw_freq=0 and m=0
      ! FOR FILTERING: exclude those waves that were damped by diffussion and for which we set gw_freq=0 and m=0
      IF (gw_frq(i,l,k) .ne. 0._r8) then           
	lambda_z(i,l, k)= (2._r8*pi)/m(i,l,k)

        if (present(kwvrdg)) then
         lambda_h_rdg(i)= (2._r8*pi)/kwvrdg(i)
         lambda_ratio(i,l,k)=lambda_z(i,l,k)/lambda_h_rdg(i)
        else
         lambda_ratio(i,l,k)=lambda_z(i,l,k)/lambda_h
	endif

	!use MF (m2/s2) to compute T' (K) using the polar eqs and dispersion 
        !rels for mid-freq Gws (see e.g. see Ern et al 2004)
        g_NT_sq(i,k)= ( gravit/(ni(i,k)*ti(i,k)) )**2.
        gw_t(i,l,k)= ( mom_flux(i,l,k)/(0.5*lambda_ratio(i,l,k)*g_NT_sq(i,k)) )**0.5 ! MF and T' are computed at interfaces (k+1/2)

        !compute Var(T')
        var_t(i,k)= var_t(i,k)+ 0.5*gw_t(i,l,k)**2.
     ELSE !set contribution from waves with gw_freq=0 to total variance zero
        var_t(i,k)= var_t(i,k)+ 0._r8
     ENDIF
  enddo 
 enddo 
enddo 

 !Compute Dt/Dz (K/m)
   do k = pver-1,1,-1
      dtdz(:,k)=(t(:,k)-t(:,k+1))/(zm(:,k)-zm(:,k+1)) !we are using t at mid-points so dtdz is computed 
      				    	 	      !at interfaces where Var(T') is defined 
   enddo

 !Compute K_eff (total effective diffusvity due to waves + Kzz)
  do k=2,pver !here we are at interfaces
     lapse_rate_sq(:,k)= (gamma_ad+dtdz(:,k-1))**2. 

     f_n_gammad(:,k)=( 1-(4._r8/3._r8)*(coriolis_f/ni(:,k))**0.5 )*gamma_ad
     b(:,k)=(cp_r-1._r8)*(var_t(:,k)/lapse_rate_sq(:,k))*( (2._r8*f_n_gammad(:,k)*(coriolis_f)**0.5 )/ti(:,k) )
     b_sq(:,k)=b(:,k)**2.
     f_ln_nf(:,k)=coriolis_f*log(ni(:,k)/coriolis_f)
     c(:,k)= (var_t(:,k)*f_ln_nf(:,k)) / (2._r8*lapse_rate_sq(:,k))

     k_eff_sqr(:,k)=b(:,k)+(b_sq(:,k)+c(:,k)+egwdffi(:,k))**0.5
     k_eff(:,k)=(k_eff_sqr(:,k))**2.
 enddo



end subroutine compute_VarT_Keff
end module gw_chem
