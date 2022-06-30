module micro_prm

implicit none

! for 30 bin aerosol and 40 bin haze/cloud, mass ratio of 2 and 2

integer, parameter :: naerosol=33, ncloud = 40
! change bin number also need to change microphysics

integer, parameter :: maindex=1, mcindex = 1 ! mass ratio = 2^(1/index)
real, parameter :: mrvalue = 1.4142135623 ! mass ratio for debug usage

integer, parameter :: nmicro_fields = 1 + naerosol + ncloud*2

real, parameter :: raer0 = 1.0e-4 ! [cm] largest droplet size
real, parameter :: rcloud0 = 1.0e-5 ! [cm] smallest droplet size

real :: Na0, ra0, siga0

real :: wall_loss_time_scale

real :: m1dropnew

logical :: docoll

integer :: distype, aer_type

logical :: inject_aer

integer :: inject_aer_step

real :: aer_source

integer :: kcdrop ! cutoff for haze and cloud droplets

integer :: krdrop ! cutoff for cloud droplets and rain

real, parameter :: rcdrop = 1.0e-4 ! cm min cloud droplet radius

real, parameter :: rrdrop = 5.0e-2 ! cm max cloud droplet radius

real, parameter :: rho_air = 1.225e-3 ! g/cm3 air density

real, parameter :: rho_H2O = 1.0 ! g/cm3 water density

real, parameter :: mol_H2O = 18.0 ! molecular weight

!real, parameter :: Rv = 461.5 ! J kg-1 K-1 already exist

real, parameter :: T0 = 293.15 ! K Temperature

real, parameter :: P0 = 1013.25 ! mb Pressure

real, parameter :: G0 = 980.0 ! cm/s2 gravity accelerator

real, parameter :: sigma_H20 = 72.8 ! dynes/cm water surface tension

integer, parameter :: iceprocs = 0 ! if iceproces=1 it is ice microphysics

! parameters to calculate critical supersaturation

real, parameter :: eps = 2.0e-16 ! for comparison
real, parameter :: epsL = 1.0e-12
real, parameter :: epsR = 1.0e-6
real, parameter :: epsS = 1.0e-38 ! for collision only
real, parameter :: eps60 = 1.0e-60 ! for zero judge
real, dimension(2), parameter :: i_vant_Hoff = (/2.0, 3.0/)
real, dimension(2), parameter :: mol_weight = (/58.443, 132.14/)
real, dimension(2), parameter :: rho_sol = (/2.17, 1.77/)
real, dimension(2), parameter :: aer_delRH = (/0.75, 0.80/)

! parameters to calculate collision
real, parameter :: eta0 = 0.0001818 ! g/cm/s
real, parameter :: delta_rho = 1e0 ! g/cm3
real, parameter :: L0 = 6.62e-6 ! cm

logical :: docs ! consider curvature and solute effect
logical :: domono ! mono-disperse aerosol
integer :: ka_mono ! mono-disperse aerosol index
real    :: ra_mono ! mono-disperse aerosol radius
real    :: ma_mono ! mono-disperse aerosol mass
real, parameter :: r0_cor = 1.86e-4 ! kinetic growth correction
integer, parameter :: DEBUG = 0
end module micro_prm
