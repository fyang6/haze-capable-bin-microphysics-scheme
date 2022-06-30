module micro_prm

implicit none

! for 30 bin aerosol and 60 bin cloud, mass ratio of 2 and (2)^(0.5), the range
! is from about 1 nm to 1 mm.

integer, parameter :: naerosol=33, ncloud = 33
! change bin number also need to change microphysics

integer, parameter :: maindex=1, mcindex = 1 ! mass ratio = 2^(1/index)
real, parameter :: mrvalue = 1.4142135623 ! mass ratio for debug usage

integer, parameter :: nmicro_fields = 1 + naerosol + ncloud*2

real, parameter :: rcloud0 = 1.0e-4 ! [cm] smallest droplet size

real :: Na0, ra0, siga0

real :: wall_loss_time_scale

real :: m1dropnew

logical :: docoll

integer :: distype, aer_type

logical :: inject_aer

integer :: inject_aer_step

real, parameter :: rho_air = 1.225e-3 ! g/cm3 air density

real, parameter :: rho_H2O = 1.0 ! g/cm3 water density

real, parameter :: mol_H2O = 18.02 ! molecular weight

!real, parameter :: Rv = 461.5 ! J kg-1 K-1 already exist

real, parameter :: T0 = 293.15 ! K Temperature

real, parameter :: P0 = 1013.25 ! mb Pressure

real, parameter :: G0 = 980.0 ! cm/s2 gravity accelerator

real, parameter :: sigma_H20 = 72.8 ! dynes/cm water surface tension

integer, parameter :: iceprocs = 0 ! if iceproces=1 it is ice microphysics

integer, parameter :: krdrop = 20 ! cutoff for cloud droplets and rain

! parameters to calculate critical supersaturation

real, parameter :: eps = 2.0e-16 ! for comparison
real, parameter :: epsL = 1.0e-12
real, parameter :: epsR = 1.0e-6
!real, parameter :: epsS = 1.0e-38 ! tiny tiny threshold
real, parameter :: epsS = 1.0e-38 ! for collision only
real, parameter :: eps60 = 1.0e-60 ! for zero judge
integer, dimension(2), parameter :: i_vant_Hoff = (/2, 3/)
real, dimension(2), parameter :: mol_weight = (/58.44, 132.14/)
real, dimension(2), parameter :: rho_sol = (/2.16, 1.77/)


! parameters to calculate collision
real, parameter :: eta0 = 0.0001818 ! g/cm/s
real, parameter :: delta_rho = 1e0 ! g/cm3
real, parameter :: L0 = 6.62e-6 ! cm
end module micro_prm
