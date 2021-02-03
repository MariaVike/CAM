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
           t, rhoi, nm, ni, c, tau, egwdffi, ubi, k_wave, xi, gw_enflux,  &
           zm, zi)
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
  ! Pressure coordinates.
  type(Coords1D), intent(in) :: p
  ! Time step.
  real(r8), intent(in) :: dt
  ! Midpoint temperature.
  real(r8), intent(in) :: t(ncol,pver)
  ! Interface density (Kg m-3)
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at interfaces.
  real(r8), intent(in) :: ubi(ncol,pver+1)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Wave Reynolds stress.
  real(r8), intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(in) :: egwdffi(ncol,pver+1)
  ! Midpoint and Interface altitudes above ground (m).
  real(r8), intent(in) :: zm(ncol,pver), zi(ncol,pver+1)    


  real(r8), intent(out) :: k_wave(ncol,pver) !total over entire wave spectrum for each GW source (i.e. Beres and C&M)
  real(r8), intent(out) :: xi(ncol,pver)
  real(r8), intent(out) :: gw_enflux(ncol,pver)

  
  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  integer :: k, l, i
  ! Gas consant for dry air (m2 K-1 s-2)
  real(r8), parameter :: R_air = 287._r8
  ! Adiabatic lapse rate and R/cp ratio
  real(r8) :: gamma_ad, r_cp 
  ! Interface temperature.
  real(r8) :: ti(ncol,pver+1)
  ! The vertical wavelength and the lmbd_h/lmbd_z ratio
  real(r8) :: lambda_z(ncol,-band%ngwv:band%ngwv,pver+1) 
  real(r8) :: lambda_ratio(ncol,-band%ngwv:band%ngwv,pver+1) 
  ! The absolute momenum flux compute from tau
  real(r8) :: mom_flux(ncol,-band%ngwv:band%ngwv,pver+1)
  ! GW intrinsic frequency 
  real(r8):: gw_frq(ncol,-band%ngwv:band%ngwv, pver+1)
  ! GW intrinsic phase speed
  real(r8):: c_i(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Vertical wavenumber m
  real(r8):: m(ncol,-band%ngwv:band%ngwv,pver+1) 
  ! GW temperature perturbation
  real(r8):: gw_t(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Variance of dT'/dz
  real(r8):: var_t(ncol,pver+1)
  ! Pressure scale height
  real(r8):: Hp(ncol,pver+1)
  ! energy term part of the computation for the energy flux
  real(r8):: energy(ncol,pver+1)
  ! Temporary values used for calculations
  real(r8), dimension(ncol, -band%ngwv:band%ngwv,pver+1) ::  frq_n, frq_m, &
				      			     m_sq, gw_t_sq, a
  real(r8), dimension(ncol,pver+1) :: g_NT_sq, b
  real(r8), dimension(ncol, pver)  :: dtdp, dtdz, one_plus_xi, one_min_xi

 logical  :: pressure_coords = .False. 

 !compute adiabatic lapserate and R/Cp ratio
 gamma_ad=gravit/cpair
 r_cp= R_air/cpair
 !compute temperature at interface (call function)
  ti(:,2:pver)=midpoint_interp(t)
 
 !compute a few quantities needed for the computation of K_wave
do i=1,ncol
 do k = 1, pver+1 
  do l = -band%ngwv, band%ngwv
        !compute the momentum flux MF (m2/s2) from the wave stress tau (Pa)
         mom_flux(i,l,k)=tau(i,l,k)/rhoi(i,k)     

	 !compute gw intrinsic frequency and vertical wavenumber
         c_i(i,l,k)=c(i,l)-ubi(i,k) 

         if (c_i(i,l,k) .ne. 0) then
           gw_frq(i,l,k)=abs(c_i(i,l,k))/lambda_h
           m(i,l,k)= (ni(i,k)*band%kwv)/gw_frq(i,l,k)
         else
         ! at critical levels where c_i=0 (i.e. c=ubi) tau=0
	 ! set everything else to zero too otherwise in the
	 ! computation of m we get division by zeros and generation of NaNs and Inf
          gw_frq(i,l,k)=0._r8  
	  m(i,l,k)=0._r8
         endif
      
  enddo
 enddo
enddo


 !compute Xi_inst = Var(dT'/dz)/(gamma_ad + dT/dz)**2
do i=1,ncol
 do k = 1, pver+1
  var_t(i,k)=0._r8
  do l = -band%ngwv, band%ngwv

      !as before for critical levels where c_i=0 and thus gw_freq=0 
       if (gw_frq(i,l,k) .ne. 0._r8) then           
	lambda_z(i,l, k)= (2._r8*pi)/m(i,l,k)
        lambda_ratio(i,l, k)=lambda_z(i,l, k)/lambda_h

	!use MF (m2/s2) to compute T' (K) using the polar eqs and dispersion 
        !rels for mid-freq Gws (see e.g. see Ern et al 2004)
        g_NT_sq(i,k)= ( gravit/(ni(i,k)*ti(i,k)) )**2.
        gw_t(i,l,k)= ( mom_flux(i,l,k)/(0.5*lambda_ratio(i,l,k)*g_NT_sq(i,k)) )**0.5 ! MF and T' are computed at interfaces (k+1/2)

        !compute Var(dT'/dz)
        var_t(i,k)= var_t(i,k)+ (m(i,l,k)**2.)*(gw_t(i,l,k)**2.)*0.5
     else
        var_t(i,k)= var_t(i,k)+ 0._r8
     endif
           

	!if (masterproc) then 
        !if (l .eq. 0) then
	!write(iulog,*)  'i,l,k =', i,l,k
        !write(iulog,*) 'TAU in gw_chem', tau(:,l,k)
        !write(iulog,*) 'MF', mom_flux(:,l,k)
        !write(iulog,*) 'gw_t', gw_t(:,l,k)
        !write(iulog,*) 'c', c(:,l)
        !write(iulog,*) 'ubi', ubi(:,k)
        !write(iulog,*) 'ci', c_i(:,l,k)
        !write(iulog,*) 'm', m(:,l,k)
	!write(iulog,*) 'band%kwv', band%kwv
	!write(iulog,*)  'var_t', var_t(:,k)
        !!write(iulog,*)  'm2', (m(:,l,k)**2)
        !write(iulog,*)  '(gw_t2)*0.5',(gw_t(:,l,k)**2)*0.5 
        !write(iulog,*)  'var_t + m2*(gw_t2)*0.5', var_t(:,k)+ (m(:,l,k)**2)*(gw_t(:,l,k)**2)*0.5  
        !endif
	!endif

  enddo 
 enddo 
enddo 

 !use model pressure coords
 if (pressure_coords) then
   do k = 2, pver
      dtdp(:,k)=(t(:,k)-t(:,k-1)) * p%rdst(:,k-1)    !using model pressure coords (p%rdst=1/Delta_p)
      xi(:,k)= var_t(:,k)/(gamma_ad+dtdp(:,k))**2.  !we are using t at mid-points so value is at interface
   enddo
 else
 !or z coordinates
   do k = pver-1,1,-1
      dtdz(:,k)=(t(:,k)-t(:,k+1))/(zm(:,k)-zm(:,k+1))
      xi(:,k)= var_t(:,k)/(gamma_ad+dtdz(:,k))**2.
   enddo
 endif

 !compute energy_flux
do i=1,ncol
 do k = 1, pver+1
  energy(i,k)=0._r8
  Hp(i,k)= (R_air*ti(i,k))/gravit
  do l = -band%ngwv, band%ngwv
     !as before for c_i=0
     if (gw_frq(i,l,k) .ne. 0._r8) then   
       frq_n(i,l,k)=gw_frq(i,l,k)**2./ni(i,k)**2.
       frq_m(i,l,k)=gw_frq(i,l,k)/m(i,l,k)    
       m_sq(i,l,k)=m(i,l,k)**2.
       gw_t_sq(i,l,k)=gw_t(i,l,k)**2.
       a(i,l,k)=(1._r8-2._r8*r_cp*(gw_frq(i,l,k))**2.)/ni(i,k)**2.
       b(i,k)=2._r8*Hp(i,k)
       energy(i,k)=energy(i,k) + (1-frq_n(i,l,k))*frq_m(i,l,k)* & 
	        ( (m_sq(i,l,k)*gw_t_sq(i,l,k)*0.5)/(m_sq(i,l,k)+ & 
	        a(i,l,k)**2./b(i,k)**2.) )
      else 
       energy(i,k)=energy(i,k) + 0._r8
      endif 
  enddo
 enddo
enddo

 if (pressure_coords) then
   do k = 2, pver  
    dtdp(:,k)=(t(:,k)-t(:,k-1)) * p%rdst(:,k-1)
    gw_enflux(:,k)= ( 1._r8/( Hp(:,k)*(gamma_ad+dtdp(:,k))**2.) ) *energy(:,k)
   enddo
 else
 !or z coordinates
   do k = pver-1,1,-1
      dtdz(:,k)=(t(:,k)-t(:,k+1))/(zm(:,k)-zm(:,k+1))
      gw_enflux(:,k)= ( 1._r8/( Hp(:,k)*(gamma_ad+dtdz(:,k))**2. ) ) *energy(:,k)
   enddo
 endif

 !Finally compute k_wave
   do k = 1, pver
      one_plus_xi(:,k)=1._r8+xi(:,k)
      one_min_xi(:,k) =1._r8-xi(:,k)    
      k_wave(:,k)=(one_plus_xi(:,k)/one_min_xi(:,k))*xi(:,k)*egwdffi(:,k)+ & 
                  ( (1._r8-r_cp)/one_min_xi(:,k) )*gw_enflux(:,k)
   enddo

   !set variables at model top (k=1) to be zero. Here tau is almost zero and
   !xi,k_wave and enflux are unphysically high (this causes a problem with the
   !netcdf output files. The writer cannot deal with such high values)
   !N.B. tau is set to zero at k=1 if tau_0_ubc=.True.
   k_wave(:,1)=0._r8
   xi(:,1)=0._r8
   gw_enflux(:,1)=0._r8

!IF (masterproc) then
!       do i=1,ncol
!        do k = 1, pver+1
!         do l = -band%ngwv, band%ngwv 
!           if (tau(i,l,k) .gt. 0) then
!		write (iulog,*)  'i,l,k =', i,l,k
!		write (iulog,*)  'ci', c_i(i,l,k)
!		write (iulog,*)  'gw_freq', gw_frq(i,l,k)
!		write (iulog,*)  'c', c(i,l)
!		write (iulog,*)  'ubi', ubi(i,k)
!          	write (iulog,*) 'TAU in gw_chem', tau(i,l,k)
!          	write (iulog,*) 'MF', mom_flux(i,l,k)
!                write (iulog,*) 'gw_t=', gw_t(i,l,k)
!          	write (iulog,*) 'var_t=', var_t(i,k)
!          	write (iulog,*) 'temp & N=',  ti(i,k), ni(i,k)
!          	write (iulog,*) 'm sq=',  m(i,l,k)
!		write (iulog,*) 'lambda_z=', lambda_z(i,l,k)
!		write (iulog,*) 'g_NT_sq=', g_NT_sq(i,k)
!          	write (iulog,*) 'dtdz', dtdz(i,k)
!		write (iulog,*) 'energy=', energy(i,k)
!		write (iulog,*) 'gw_enflux=', gw_enflux(i,k)
!		write (iulog,*) 'xi=', xi(i,k)
!		write (iulog,*) 'k_wave=', k_wave(i,k)     		
!          endif
!         enddo
!        enddo
!      enddo
!END IF 


end subroutine effective_gw_diffusivity

end module gw_chem
  
