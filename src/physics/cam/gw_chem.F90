module gw_chem

!
! This module contains code to compute the wave‐driven constituent transport 
! due to non-breaking gravity waves (Gardner et al, ESS, 2019 ; Gardner et al, JGR, 2018).
! Here we compute the effective wave diffusivity (K_wave) as a function of the 
! eddy diffusivity (Kzz) and of thevariances of the temperature and lapse rate fluctuations
!

use gw_utils,  only: r8
use coords_1d, only: Coords1D

implicit none
private
save

! Public interface.
public :: effective_gw_diffusivity

contains

!==========================================================================

subroutine effective_gw_diffusivity (ncol, band, lambda_h, p, dt,    &
           t, rhoi, nm, ni, c, tau, egwdffi, k_wave, xi, gw_enflux,  &
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
  ! Interface densities.
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
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
  integer :: k, l
  ! Gas consant for dry air (m2 K-1 s-2)
  real(r8), parameter :: R_air = 287._r8
  ! Adiabatic lapse rate 
  real(r8), parameter :: gamma_ad=gravit/cpair
  real(r8), parameter :: r_cp= R_air/cpair
  ! Interface temperature.
  real(r8) :: ti(ncol,pver+1)
  ! The vertical wavelength and the lmbd_h/lmbd_z ratio
  real(r8) :: lambda_z, lambda_ratio
  ! The absolute momenum flux compute from tau
  real(r8) :: MF(ncol,-band%ngwv:band%ngwv,pver+1)
  ! GW intrinsic frequency 
  real(r8):: gw_frq(ncol,-band%ngwv:band%ngwv) 
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
  ! Temporary values for computing vertical derviatives
  real(r8):: dtdp, dtdz
  ! Temporary values used for calculations
  real(r8) :: frq_n, frq_m, m_sq, gw_t_sq 
  real(r8) :: a, b, one_plus_xi, one_min_xi
     

 logical  :: pressure_coords = .True. 

 !compute temperature at interface (call function)
  ti(:,2:pver)=midpoint_interp(t)
 
 !compute a few quantities needed for the computation of K_wave
 do k = 1, pver+1 
  do l = -band%ngwv, band%ngwv
        !compute the momentum flux MF from the wave stress tau
        MF(:,l,k)=tau(:,l,k)/rhoi(:,k)

        gw_frq(:,l)=c(:,l)/lambda_h
        m(:,l,k)=(ni(:,k)**2 + band%kwv**2 /gw_frq(:,l)**2)**0.5        
  enddo
 enddo

 !compute Xi_inst = Var(dT'/dz)/(gamma_ad + dT/dz)**2
 do k = 1, pver+1
  var_t=0.
  do l = -band%ngwv, band%ngwv
        !use MF to compute T' using the polar eqs and dispersion 
        !rels for mid-freq Gws (see e.g. see Ern et al 2004)
        lambda_z= 2._r8*pi/m(:,l, k)
        lambda_ratio=lambda_z/lambda_h
        gw_t(:,l,k)= MF(:,l,k)/(0.5*rhoi(:,k)*lambda_ratio*(gravit/ni(:,k)*ti(:,k)) )**0.5 ! MF and T' are computed at interfaces (k+1/2)
        !compute Var(dT'/dz)
        var_t(:,k)= var_t(:,k)+ (m(:,l,k)**2)*(gw_t(:,l,k)**2)*0.5        
  enddo 
 enddo      


 !use model pressure coords
 if (pressure_coords) then
   do k = 2, pver
      dtdp=t(:,k)-t(:,k-1) * p%rdst(:,k-1)    !using model pressure coords (p%rdst=1/Delta_p)
      xi(:,k)= var_t(:,k)/(gamma_ad+dtdp)**2  !we are using t at mid-points so value is at interface
   enddo
 else
 !or z coordinates
   do k = pver,1,-1
      dtdz=t(:,k)-t(:,k+1)/zm(:,k)-zm(:,k+1)
      xi(:,k)= var_t(:,k)/(gamma_ad+dtdz)**2
   enddo
 endif

 !compute energy_flux
 do k = 1, pver+1
  energy=0.
  Hp(:,k)= R_air*ti(:,k)/gravit
  do l = -band%ngwv, band%ngwv
     frq_n=gw_freq(:,l)**2/ni(:,k)**2
     frq_m=gw_freq(:,l)/m(:,l,k)
     m_sq=m(:,l,k)**2
     gw_t_sq=gw_t(:,l,k)**2
     a=(1-2*r_cp*gw_freq(:,l)**2)/ni(:,k)**2
     b=2*Hp(:,k)
     energy(:,k)=energy(:,k) + (1-frq_n)*frq_m*( (m_sq*gw_t_sq*0.5)/(m_sq+a**2/b**2) )
  enddo
 enddo

 if (pressure_coords) then
   do k = 2, pver  
    dtdp=t(:,k)-t(:,k-1) * p%rdst(:,k-1)
    gw_enflux(:,k)= (1/Hp(:,k)*gamma_ad+dtdp**2)*energy(:,k)
   enddo
 else
 !or z coordinates
   do k = pver,1,-1
      dtdz=t(:,k)-t(:,k+1)/zm(:,k)-zm(:,k+1)
      gw_enflux(:,k)= (1/Hp(:,k)*gamma_ad+dtdz**2)*energy(:,k)
   enddo
 endif

 !Finally compute k_wave
 do k = 1, pver
    one_plus_xi=1+xi(:,k)
    one_min_xi =1-xi(:,k)    
    k_wave(:,k)=(one_plus_xi/one_min_xi)*xi(:,k)*egwdffi(:,k)+ & 
                (1-r_cp/one_min_xi)*gw_enflux(:,k)
 enddo
  

end subroutine effective_gw_diffusivity

end module gw_chem
  
