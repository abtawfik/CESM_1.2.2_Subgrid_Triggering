module hcf_trigger

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: hcf_trigger
!
! !DESCRIPTION:
! Call a new test for convection based on the buoyant condensation level (BCL)
! This method assumes moist deep convection is triggered when saturation occurs
! at the top of the mixed layer.  Mixed layer is defined by the neutral buoyancy
! level of a surface parcel (usually 2-meter parcels).
! Reference: Tawfik and Dirmeyer 2013 GRL: A processed based framework for
!            quantifying the atmospheric background state of surface triggered
!            convection

!
! !USES:
  use shr_kind_mod,    only : r8 => shr_kind_r8
  use spmd_utils  ,    only : masterproc
  use ppgrid      ,    only : pcols, pver, pverp
  use spmd_utils  ,    only : masterproc
!
! !PUBLIC TYPES:
  implicit none

! !PUBLIC MEMBER FUNCTIONS:
!
  public :: hcf_conv_trigger           ! Interface that preps input variables for calc conv threshold
  public :: mixed_layer                ! computes threshold for a single column

!------------------------------------------------------------------------------

   contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  hcf_conv_trigger
!
!
!
! !INTERFACE:
  subroutine hcf_conv_trigger( ncol, tmp_k, press, qhum, tbm )
!
!
!
! !DESCRIPTION:
! Call a new test for convection based on the buoyant condensation level (BCL)
! This method assumes moist deep convection is triggered when saturation occurs
! at the top of the mixed layer.  Mixed layer is defined by the neutral buoyancy
! level of a surface parcel (usually 2-meter parcels).
! Reference: Tawfik and Dirmeyer 2013 GRL: A processed based framework for
!            quantifying the atmospheric background state of surface triggered
!            convection
!**** Input variables ****
!--- press  =  pressure levels in [Pascals]
!--- temp   =  potential temperature profile in [Kelvin]
!--- tmp_k  =  temperature profile in [Kelvin]
!--- qhum   =  specific humidity profile in [kg H2O / kg air]
!--- nnn    =  number of vertical levels
!
!**** Outut variables ****
!--- tbm    =  potential temperature of convective threshold [K]
!
   integer ,                        intent(in ) :: ncol
   real(r8), dimension(pcols,pver), intent(in ) :: qhum, tmp_k, press
   real(r8), dimension(pcols)     , intent(out) :: tbm

   real(r8), dimension(pcols,pver) :: temp, pres
   real(r8), dimension(pcols,pver) :: temp1, pres1, tmpk1, qhum1
   real(r8), dimension(pver)       :: temp0, pres0, tmpk0, qhum0
   integer                         :: ii, k, nnn, aa, bb
!------------------------------------------------------------------------------


   !*** Calculate potential temperature for all layers
   pres(:ncol,:)  =  press(:ncol,:) * 100._r8         !-- convert to Pa from hPa


   !*** Rearrange to bottom to top
   ii   = 1
   nnn  = pver-1     !--- number of levels
   do k = pver,1,-1
      temp1(:ncol,ii)  =  temp (:ncol,k)
      tmpk1(:ncol,ii)  =  tmp_k(:ncol,k)
      qhum1(:ncol,ii)  =  qhum (:ncol,k)
      pres1(:ncol,ii)  =  pres (:ncol,k)
      ii = ii + 1
   end do



   !*******************************************
   !*** Calculate heated condensation level ***
   !*******************************************
   do ii = 1,ncol

          temp0(:nnn)  =  temp1(ii,:nnn)
          tmpk0(:nnn)  =  tmpk1(ii,:nnn)
          pres0(:nnn)  =  pres1(ii,:nnn)
          qhum0(:nnn)  =  qhum1(ii,:nnn)

          call mixed_layer ( nnn, tmpk0(:nnn), pres0(:nnn), qhum0(:nnn), tbm(ii), .false.)

   end do


!------------------------------------------------------------------------------
end subroutine hcf_conv_trigger
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  mixed_layer
!
!
!
! !INTERFACE:
  subroutine mixed_layer( nlev1, tmp_in, press_in, qhum_in, TBM, BUGGIN )
!
!
! !DESCRIPTION
! Call a new test for convection based on the buoyant condensation level (BCL)
! This method assumes moist deep convection is triggered when saturation occurs
! at the top of the mixed layer.  Mixed layer is defined by the neutral buoyancy 
! level of a surface parcel gradually mixed due to incremental surface heating.
!
! !NOTES
! Vertical profiles are assumed to be bottom-up (surface to tropopause)
! Assumes first model level is either 2-m or surface layer
!
! !REFRENCE 
! ***If you use this routine please reference the publication below***
! Tawfik and Dirmeyer 2014 GRL: A processed based framework for 
! quantifying the atmospheric background state of surface triggered
! convection
!
! !Input variables
! press  =  pressure levels in [Pascals]
! tmp_k  =  temperature profile in [Kelvin]
! qhum   =  specific humidity profile in [kg H2O / kg air]
! n      =  number of vertical levels 
!
! !Output variables
! tbm  =  potential temperature threshold for convective initiation
!                                                                           
!                                                                           
   implicit none
!  
! Input/Output Variables
!  
   integer , intent(in   )                       ::  nlev1
   real(r8), intent(in   ), dimension(nlev1)     ::  tmp_in, qhum_in, press_in
   real(r8), intent(out  )                       ::  TBM
   logical , intent(in   )                       ::  BUGGIN
!  
! Local variables  
!    
   real(r8), parameter   ::  p_ref=1e5 , Lv=2.5e6, cp=1005.7, R_cp=287.04/1005.7
   real(r8), parameter   ::  grav  = 9.81, Rd=287.04
   real(r8), parameter   ::  by100 = 100. 
   real(r8), parameter   ::  large_positive = 99999.       ! used to replace missing value when min bound is needed   
   real(r8), parameter   ::  large_negative = -99999.      ! used to replace missing value when max bound is needed 
   real(r8), parameter   ::  t0=273.15, ep=0.622, es0=6.11, a=17.269, b=35.86
   integer               ::  nlatlon
   integer               ::  zz
   integer               ::  nlev

   logical, dimension(nlev1)    ::  notmissing
   integer, dimension(nlev1)    ::  idef
   real(r8), dimension(nlev1)   ::  rhoh

   real(r8), dimension(nlev1)   ::  allmissing
   real(r8), dimension(nlev1)   ::  pbar, qdef, qmix, qsat, dpress, qbar, logp, tbar
   real(r8), dimension(nlev1)   ::  tmp_k, press
   real(r8), dimension(nlev1)   ::  qhum

   real(r8)                     ::  miny, maxy
   real(r8), dimension(nlev1)   ::  min0, max0, dummy
   real(r8), dimension(nlev1)   ::  p_up0, t_up0, h_up0, q_up0, p_lo0, t_lo0, h_lo0, q_lo0
   real(r8)                     ::  p_up, t_up, h_up, q_up, p_lo, t_lo, h_lo, q_lo
   real(r8), dimension(nlev1)   ::  dufus 
   real(r8)                     ::  onemep
   real(r8)                     ::  BCLP
   logical                      ::  transition

!-----------------------------------------------------------------------------  

      !-- Store temporary working arrays and initialize 
      nlev   =  nlev1
      tmp_k  =  tmp_in
      qhum   =  qhum_in
      press  =  press_in


      !-- Initialize middle level variables 
      notmissing   =  .false.

      !-- Get Level indices 
      do zz = 1,nlev
         idef(zz)  =  zz
      end do


      !-- Calculate middle layer temperature 
      tbar(1:nlev-1)    =  ( tmp_k(2:nlev  )*log(press(2:nlev  ))   +  &
                             tmp_k(1:nlev-1)*log(press(1:nlev-1)) ) /  &
                             log(press(2:nlev) * press(1:nlev-1))

      !-- Calculate pressure difference between layers     
      dpress(1:nlev-1)  =  press(1:nlev-1) - press(2:nlev)

      !-- Calculate middle layer specific humidity (kg/kg)  
      qbar(1:nlev-1)    =  ((qhum(2:nlev  )*log(press(2:nlev  ))  + &
                             qhum(1:nlev-1)*log(press(1:nlev-1))) / &
                             log(press(2:nlev)* press(1:nlev-1)))

      !-- Calculate middle layer pressure (Pa)   
      pbar(1:nlev-1)    =  press(1:nlev-1)

      !-- Calculate log of pressure     
      logp(1:nlev-1)    =  log(pbar(1:nlev-1))

      !-- Calculate calculate mixed moisture (kg/m2) and column density (kg/m2) 
      qmix(1:nlev-1)  =  qbar  (1:nlev-1)*dpress(1:nlev-1)/grav
      rhoh(1:nlev-1)  =  dpress(1:nlev-1)/grav

      !-- Calculate calculate mixed humidity (kg/kg)  
      do zz = 2,nlev-1
         qmix(zz)  =  qmix(zz-1) + qmix(zz)
         rhoh(zz)  =  rhoh(zz-1) + rhoh(zz)
      end do

      !-- Calucalte the saturation specific humidity (kg/kg)  
      pbar(1:nlev-1)  =  pbar(1:nlev-1)/1e2
      onemep =  1.0 - ep
      qsat(1:nlev-1)  =  by100*0.01*(ep* (es0*exp((a*( tbar(1:nlev-1)-t0))/( tbar(1:nlev-1)-b))) ) /  &
                         (pbar(1:nlev-1)-onemep*(es0*exp((a*( tbar(1:nlev-1)-t0))/( tbar(1:nlev-1)-b))))
      qsat(1:nlev-1)  =  qsat(1:nlev-1)/(1.+qsat(1:nlev-1))
      pbar(1:nlev-1)  =  pbar(1:nlev-1)*1e2

      qmix(1:nlev-1)  =  qmix(1:nlev-1) / rhoh(1:nlev-1)
      qdef(1:nlev-1)  =  qsat(1:nlev-1) - qmix(1:nlev-1)

 

      !--- Check to make sure qdef is always negative when outside of the tropo 
      !--- Assume a tropopause height of 10 km; so BCL cannot be higher  
      where( pbar.le.5000. )
         qdef  =  -1.
      endwhere



      !***********************************************************  
      !***   Calculate slope of each variable to find the      ***  
      !***   y-intercept for each variable;                    ***  
      !***   Meaning locate the two data points surrounding    ***  
      !***   the sign change in qdef and linearly interpolate  ***  
      !***   to find the "zero point"                          ***  
      !***********************************************************                 
      !----- Get the higher (from ground) level point in the qdef sign transition  
      !----- Find the point where the sign first turns negative from the ground up 
      transition  =  .false.
      do zz = 2,nlev-1
         if( qdef(zz-1).gt.0 .and. qdef(zz).le.0 ) then
            transition  =  .true.
            p_up  =  logp(zz)
            t_up  =  tbar(zz)
            q_up  =  qdef(zz)

            p_lo  =  logp(zz-1)
            t_lo  =  tbar(zz-1)
            q_lo  =  qdef(zz-1)
            exit
         end if
      end do

      !--- If there is no BCL height (e.g. no transition in saturation such as can occur 
      !--- over the Antarctic) then set values to heightest level
      if(.not.transition) then
         BCLP   =  exp( logp(nlev-1) )
         TBM    =  tbar(nlev-1) * ((p_ref/BCLP)**(R_cp))
         return
      end if

      !--- Calculate output variables; BCL height, BCL pressure,
      !--- Buoyant Mixing Potential Temp, and Potential Temperature Deficit 
      !*** Note ***  
      !--- !!!!  This is the Virtual Potential Temperature (K) !!!  
      BCLP   =  exp( p_up - ((p_up-p_lo)/(q_up-q_lo))*q_up )
      TBM    =     ( t_up - ((t_up-t_lo)/(q_up-q_lo))*q_up ) * ((p_ref/BCLP)**(R_cp))


end subroutine mixed_layer



end module hcf_trigger


