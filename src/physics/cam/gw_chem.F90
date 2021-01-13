module gw_chem

!
! This module contains code to compute the transport of 
! chemical species due to gravity wave propagation (Gardner et al, JGR, 2018)
!

use gw_utils,  only: r8
use coords_1d, only: Coords1D

implicit none
private
save

! Public interface.
public :: effective_gw_diffusivity

! gas consant for dry air (m2 K-1 s-2)
real(r8), parameter :: R_air = 287._r8

contains

!==========================================================================

subroutine effective_gw_diffusivity (ncol, band, lambda_h, p, dt,      &
              t, rhoi, nm, ni, c, tau, egwdffi, k_wave, xi, gw_enflux, &
              k_wave_atc, xi_atc, gw_enflux_atc, egwdffi_atc)
!-----------------------------------------------------------------------
! Compute K_wave (wave effective diffusivity) etc...
! ....
!-----------------------------------------------------------------------
use gw_common, only: GWBand, pver
use physconst, only: cpair, cpairv, gravit
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
  ! Midpoint and interface temperatures.
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


  real(r8), intent(out) :: k_wave(ncol,pver+1) !total over entire wave spectrum for each GW source (i.e. Beres and C&M)
  real(r8), intent(out) :: xi(ncol,pver+1)
  real(r8), intent(out) :: gw_enflux(ncol,pver+1)

  real(r8), intent(out) :: k_wave_atc(ncol,-band%ngwv:band%ngwv,pver+1) ! functions of wave-frequency
  real(r8), intent(out) :: xi_atc(ncol,-band%ngwv:band%ngwv,pver+1)
  real(r8), intent(out) :: gw_enflux_atc(ncol,-band%ngwv:band%ngwv,pver+1)
  real(r8), intent(out) :: egwdffi_atc(ncol,-band%ngwv:band%ngwv,pver+1)

  
  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  !integer :: k, l, m, i

  
  !integer :: 

  ! 
  !real(r8) ::




  k_wave=5
  xi=1
  gw_enflux=8 

end subroutine effective_gw_diffusivity

end module gw_chem
  
