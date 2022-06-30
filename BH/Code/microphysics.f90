module microphysics

! This module is a double-moment bin microphysical schemes based on Chen and
! Lamb JAS 1994
! Reference:
! Chen, J.P. and Lamb, D., 1994. Simulation of cloud microphysical and chemical
! processes using a multicomponent framework. Part I: Description of the
! microphysical model. Journal of the Atmospheric Sciences, 51(18),
! pp.2613-2630.
!
! The code is implemented in SAM by Fan Yang @ BNL (fanyang@bnl.gov), 2021!
! History:
! version, date, description
! v0.1(F. Yang); 01/12/2021; currently only work for pi chamber simulation
! v0.2(F. Yang); 03/23/2021; adjust M after advection to avoid bin shift
! v0.3(F. Yang); 06/21/2021; bug fixed
! v0.4(F. Yang); 07/03/2021; this version includes haze particles
! v0.5(F. Yang); 07/14/2021; add statistic output for haze and activation rate

use grid, only: nx, ny, nzm, nz, masterproc, nsubdomains, RUN3D, &
              & dimx1_s, dimx2_s, dimy1_s, dimy2_s

use params, only: cp, ggr, rgas, rv, lsub, muelq, pi

use params, only: doprecip, docloud

use micro_prm

use module_ntubm

implicit none

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

! flag for 33 aerosol number and 40 cloud mass and number bins
! Adjust nmicro_fields
integer, parameter :: flag_wmass(nmicro_fields) = (/1, &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &  ! aerosol number
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &  ! cloud mass
	 0,0,0,0,0,0,0, &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
	 0,0,0,0,0,0,0 /)   ! cloud number

integer, parameter :: flag_precip(nmicro_fields) = (/0, &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &  ! aerosol
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &  ! cloud mass
	 0,0,0,0,0,0,0, &
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
	 0,0,0,0,0,0,0 /)   ! cloud number
      
integer, parameter :: index_water_vapor = 1 ! index for variables that contains water vapor
integer, parameter :: index_cloud_ice = -1 ! index for cloud ice (sedimentation)

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface fluxes
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary fluxes 
real fluxlmk (ny, nzm, 1:nmicro_fields)
real fluxrmk (ny, nzm, 1:nmicro_fields)
real fluxqmk (nx, nzm, 1:nmicro_fields)
real fluxhmk (nx, nzm, 1:nmicro_fields)

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mkdiff(nz,1:nmicro_fields) ! tendency due to vertical diffusion
real mklsadv(nz,1:nmicro_fields) ! tendency due to large-scale vertical advection

! prognostic variables:
real qt(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm) ! [g/kg] total water mixing ratio
real Naer(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,naerosol) ! [#/kg/bin] aerosol number concentration in each bin
real Ncld(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncloud) ! [#/kg/bin] cloud droplet number concentration in each bin
real Mcld(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,ncloud) ! [kg/kg/bin] cloud droplet mass in each bin 

! statistics
real qh(nx,ny,nzm) ! haze water mixing ratio [kg/kg]
real qc(nx,ny,nzm) ! liquid mixing ratio [kg/kg]
real qr(nx,ny,nzm) ! rain/drizzle mixing ratio [kg/kg]
real qlwc(nx,ny,nzm) ! total condensed water mixing ratio [kg/kg]
real qna(nx,ny,nzm) ! CCN number concentration [/mg]
real qnh(nx,ny,nzm) ! haze number concentration [/mg]
real qnc(nx,ny,nzm) ! droplet number concentration [/mg]
real qnr(nx,ny,nzm) ! rain/drizzle number concentration [/mg]
real qnt(nx,ny,nzm) ! total number concentration [/mg]

real ssatw(nx,ny,nzm) ! supersaturation over liquid water [%]

real reffc(nx,ny,nzm)
real reffi(nx,ny,nzm)
real rmean(nx,ny,nzm)

real act_rate(nx,ny,nzm)      ! activation rate [#/mg/s]
real deact_rate(nx,ny,nzm)    ! deactivation rate [#/mg/s]

real rcmean(nx,ny,nzm)    ! mean cloud droplet radius [um]
real rcsig(nx,ny,nzm)     ! standard deviation of cloud droplet radius [um]
real rceff(nx,ny,nzm)     ! effective radius [um]
real rel_dis(nx,ny,nzm)   ! relative dispersion

real difful_tend(nx, ny, nzm)

! map prognostic variables onto micro_field array:
equivalence (qt(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (Naer(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,2))
equivalence (Mcld(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,naerosol+2))
equivalence (Ncld(dimx1_s,dimy1_s,1,1),micro_field(dimx1_s,dimy1_s,1,naerosol+ncloud+2))

integer ncond


!----------------------------------------------------------------------
!
          CONTAINS
!
!----------------------------------------------------------------------

subroutine micro_setparm()
  use vars
  use grid, only: masterproc
  implicit none

  integer ierr, ios, ios_missing_namelist,place_holder

  NAMELIST /MICRO_NTUBM/ &
       Na0, & ! (cm-3) Total aerosol number concentration
       ra0, & ! (cm)
       siga0, &
       docoll, &
       distype, &
       aer_type, &
       wall_loss_time_scale, &
       inject_aer, &
       inject_aer_step, &
	   aer_source, &
	   docs, &
	   domono, &
	   ka_mono

  NAMELIST /BNCUIODSBJCB/ place_holder

  mkwle = 0. 
  mkwsb = 0. 
  mkadv = 0. 
  mkdiff = 0. 
  mklsadv = 0. 
  ncond = 1

  Na0 = 50. ! default aerosol number concentration (cm-3)
  ra0 = 50.0e-3 ! default mean aerosol radius (um)
  siga0 = 1.1 ! default standard deviation
  docoll = .false. ! flag for collision coalescence
  distype = 1  ! distribution type
  aer_type = 1 ! aerosol type
  inject_aer = .false. ! flag for aerosol injection
  inject_aer_step = 0  ! initial step for aerosol injection
  aer_source = 0.25 ! aerosol injection rate (/mg/s)
  docs = .false. ! ignore curvature and solute effects
  domono = .false. ! mono-disperse aerosol
  kcdrop = 1  ! separate haze and cloud droplets
  krdrop = 20 ! separate cloud droplets and rain
  ka_mono =22 ! index for mono-disperse aerosol

  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55)
  read(55,MICRO_NTUBM,IOSTAT=ios)

  if (masterproc) then
    write(*,'(A60,F10.4)')"Initial aerosol number concentration (cm-3):", Na0
    write(*,'(A60,F10.4)'),"Initial mean aerosol radius (cm):", ra0
    write(*,'(A60,F10.4)'),"Initial standard deviation of aerosol size distribution:", siga0
    if (aer_type == 1) then
      write(*,'(A60)')"Aerosol Type: NaCl"
    elseif (aer_type == 2) then
      write(*,'(A60)')"Aerosol Type: (NH4)2SO4"
    else
      print*, "Error: Aerosol Type Not Defined"
      stop
    end if
    if (docoll) then
      write(*,'(A60)')"Turn on collision coalescence process."
    else
      write(*,'(A60)')"Turn off collison coalescence process."
    end if
	if (docs) then
      write(*,'(A60)')"Consider curvature and solute effect."
	else
      write(*,'(A60)')"Ignore curvature and solute effect."
	end if
	if (domono) then
      write(*,'(A60)')"Consider mono-disperse aerosol."
	else
      write(*,'(A60)')"Turn off mono-disperse aerosol."
	end if

  end if

  if (ios.ne.0) then
    if (ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_NTUBM namelist'
        call task_abort()
    elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '***** No MICRO_NTUBM namelist in prm file ********'
        write(*,*) '****************************************************'
    end if
  end if
end subroutine micro_setparm


subroutine micro_init()
  use grid, only: nx, ny, nzm, dt, case
  use vars
  use params
  
  implicit none

! local variable
  integer :: i,j,k,ia,ic
  real :: aratio, cratio
  integer, parameter :: ntusbm_unit = 22

  if (masterproc) then
    write(*,*) '****************************************************'
    write(*,*) '***************** INITIALIZING NTUBM ***************'
    write(*,*) '****************************************************'
  end if

! calculate category boundaries for aerosol
  aratio = 2.0**(1.0/maindex)
  cratio = 2.0**(1.0/mcindex)
  
  do ia = naerosol, 1, -1
    if (ia == naerosol) then
      raerosol(ia) = raer0*aratio**(-1.0/3.0) 
    else
      raerosol(ia) = raerosol(ia+1)*aratio**(-1.0/3.0) 
    end if
! initial number concentration
    if (distype==1) then
      numaer0(ia) = Na0/(sqrt(2*3.1416)*log(siga0)) &
              *exp(-(log(raerosol(ia)*1.0e4/ra0))**2/2.0/(log(siga0))**2) &
              * 1.0/3.0*log(aratio) 
    elseif (distype==2) then
      if (ia==22) then
        numaer0(ia) = 0.0
      else
        numaer0(ia) = 0.0
      end if
    end if
  end do

  do ic = 1, ncloud
    if (ic == 1) then
      rcloud(ic) = rcloud0
      mcloud(ic) = 4.0*pi/3.0*rho_H2O*rcloud(ic)**3
    else
      rcloud(ic) = rcloud(ic-1)*cratio**(1.0/3.0)
      mcloud(ic) = 4.0*pi/3.0*rho_H2O*rcloud(ic)**3
      mcloud2(ic-1) = mcloud(ic)
    end if
  end do
  mcloud2(ncloud) = mcloud(ncloud)*cratio

  do ic = 1, ncloud
    if (ic==1.and.rcloud(ic)>rcdrop) then
	  kcdrop = ic
	  exit
	elseif (rcloud(ic)>rcdrop) then
	  kcdrop = ic
	  exit
	end if
  end do

  do ic = kcdrop, ncloud
    if (rcloud(ic)>rrdrop) then
	  krdrop = ic - 1
	  exit
	elseif (ic==ncloud.and.rcloud(ic)<=rrdrop) then
	  if (masterproc) print*, "Warming: no rain drops available for bin setup"
	  krdrop = ncloud
	end if
  end do

! set mono-disperse aerosol size and mass
  if (domono) then
    ra_mono = raerosol(ka_mono)
	ma_mono = 4.0*pi/3.0*rho_sol(aer_type)*raerosol(ka_mono)**3
  else
    ra_mono = 0e0
	ma_mono = 0e0
  end if

  numcld0 = 0.0
  m1dropnew = 4.0*pi/3.0*rho_H2O* ((rcloud(1) + rcloud(2))/2.0)**3

  open(unit=ntusbm_unit,file="./SRC/MICRO_NTUBM/sbm_input/term_cloud.dat",form="formatted",status="old")
  read(ntusbm_unit,*) cloud_term_vel
  close(ntusbm_unit)

  open(unit=ntusbm_unit,file="./SRC/MICRO_NTUBM/sbm_input/kernel_cloud.dat",form="formatted",status="old")
  read(ntusbm_unit,*) kernel_cloud
  close(ntusbm_unit)

  if (masterproc) print*, 'Aerosol mass ratio:', aratio
  if (masterproc) print*, 'Aerosol radius (um) and number concentration (#/cm3)'
  if (masterproc) print '(3X, I3, F12.5, F12.5)', (ia, raerosol(ia)*1.0e4,numaer0(ia), ia=1,naerosol)
  if (masterproc) print*, 'Cloud mass ratio:', cratio
  if (masterproc) print*, 'Cloud droplet radius (um), Cloud droplet mass (ug), concentration (#/cm3), Terminal velocity (cm/s)'
  if (masterproc) print '(3X, I3, F12.5, F12.5, F12.5, F12.5)', (ic, rcloud(ic)*1.0e4, mcloud(ic)*1.0e6, numcld0(ic), cloud_term_vel(ic), ic=1,ncloud)
  if (masterproc.and.inject_aer) print*,'Inject aerosol at', inject_aer_step, ' steps, at a rate of', aer_source, ' #/mg/s.'

! initialization other variables
  micro_field = 0e0
  qh = 0e0
  qc = 0e0
  qr = 0e0
  qlwc = 0e0
  qcl = 0e0
  qpl = 0e0
  qci = 0e0
  qpi = 0e0
  qna = 0e0
  qnh = 0e0
  qnc = 0e0
  qnr = 0e0
  qnt = 0e0
  Ncld = 0e0
  Mcld = 0e0
  Naer = 0e0
  difful_tend = 0e0
  act_rate = 0e0
  deact_rate = 0e0
  rcmean = 0e0
  rcsig = 0e0
  rceff = 0e0
  rel_dis = 0e0

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        qt(i,j,k) = q0(k)
        tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k) + qr(i,j,k))
        qv(i,j,k) = q0(k)
      end do
    end do
  end do 

! initialization flux
  fluxbmk = 0.
  fluxtmk = 0. 
  fluxlmk = 0. 
  fluxrmk = 0. 
  fluxqmk = 0. 
  fluxhmk = 0.

  if(docloud) call micro_diagnose()

  return
end subroutine micro_init

subroutine check_M_N(if_stop,print_what)
    use vars, only: rho
    use grid, only: nstep
    character(len=*) :: print_what
    logical, INTENT(IN) :: if_stop
    integer :: i, j , k, kr, ic
    real :: cratio, M_this, N_this, M_N
  integer, parameter :: its=1, ite=nx, jts=1, &
                        jte=ny, kts=1, kte=nzm

  do k = kts, kte
    do j = jts, jte
      do i = its, ite
        do kr = 1, ncloud
            if (masterproc) then
            
                M_this = Mcld(i,j,k,kr)*rho(k)*1.0e-3
                N_this = Ncld(i,j,k,kr)*rho(k)*1.0e-6
                if (N_this > 0 ) then
                    M_N = M_this / N_this
                    if (M_N > mcloud2(kr) .or. M_N < mcloud(kr)) then
                        print*, print_what,nstep, i,j,k,kr
                        print*, 'ERROR1: M_N out of boundary'
                        if (if_stop) stop
                    end if
                else if (M_this > 0) then
                    print*, print_what,nstep, i,j,k,kr
                    print*, 'ERROR2: M not 0, M, N', M_this, N_this 
                    if (if_stop) stop
                end if 
            end if
        end do
      end do
    end do
  end do
end subroutine check_M_N

subroutine fix_MN_m1(Mrho, Nrho, MN, bd, bd2)
	real, INTENT(IN) :: Nrho, MN, bd, bd2
	real, INTENT(INOUT) :: Mrho

	real :: new_MN, bs 
	real :: m_num
	
	! modify m_num to find a good value as MN
	m_num = 0.01

	if (MN>bd) then ! right boundary
		bs = bd - bd2
		new_MN = bd - m_num * bs
	else ! left boundary
		bs = bd2 - bd
		new_MN = bd + m_num * bs
	end if
	Mrho = new_MN * Nrho
end subroutine fix_MN_m1

subroutine fix_MN_m2(Mrho, Nrho, MN, bd, bd2)
	real, INTENT(IN) :: Nrho, MN, bd, bd2
	real, INTENT(INOUT) :: Mrho

	real :: m_num, k 
	real :: x1, x2
	real :: xs ! x_star
	real :: N_temp

	m_num = 1.0e8
	k = m_num

	if (MN>bd) then ! right boundary
		x2 = bd
		x1 = x2 - sqrt(2*Nrho/k)
		if (x1 < bd2) then ! k make the whole bin covered
			x1 = bd2
			k = 2.0e0 * Nrho / (x2 - x1)**2
			if (masterproc) print*,'strong warning: fix_MN_m2 k to small',k,'new k'
		end if
		xs = x1
		N_temp = ( x2 - x1 ) * ( - k * (xs - (x1 + x2) / 2.e0))
	else ! left boundary
		x1 = bd
		x2 = sqrt(2*Nrho/k) + x1
		if (x2 > bd2) then
			x2 = bd2
			k = 2.0e0 * Nrho / (x2 - x1)**2
		end if
		xs = x2
		k = -k
		N_temp = ( x2 - x1 ) * ( - k * (xs - (x1 + x2) / 2.e0))
	end if
	call line2MN(x1, x2, k, xs, 0.0e0, Mrho, N_temp)
		
	if (abs(N_temp - Nrho)>eps) then
		print*, 'warning: fix_MN line2MN from N',Nrho,'to',N_temp,'diff',N_temp-Nrho
		print*, 'bd,bd2,x1,x2,xs,k:',bd,bd2,x1,x2,xs,k
	end if
end subroutine fix_MN_m2

subroutine fix_MN(this_M, this_N, this_rho, kr)
	real, INTENT(IN) :: this_rho
	integer, INTENT(IN) :: kr
	real, INTENT(INOUT) :: this_M, this_N

	real :: Mrho, Nrho
	real :: MN
	real :: bd, bd2
	integer :: method 
	! method 1: MN at boundary+-bin_size*m_num
	! method 2: make the |k|=m_num

	method = 1

	if (this_N>eps60) then
		Mrho = this_M * this_rho * 1.0e-3
		Nrho = this_N * this_rho * 1.0e-6
		MN = Mrho/Nrho
		if (MN > mcloud2(kr)) then
			bd = mcloud2(kr)
			bd2 = mcloud(kr)
		else if (MN < mcloud(kr)) then
			bd = mcloud(kr)
			bd2 = mcloud2(kr)
		else
			return
		end if

		if (method==1) then
			call fix_MN_m1(Mrho, Nrho, MN, bd, bd2)
		else if (method==2) then
			call fix_MN_m2(Mrho, Nrho, MN, bd, bd2)
		end if

		this_M = Mrho / this_rho / 1.0e-3
	else if (this_N<-eps60) then
		print*, 'WARNING: after advection, N=',this_N,'<0'
		print*, 'Solution: Force it to zero'
		this_N = 0.0e0
		this_M = 0.0e0
		
	else ! this_N == 0
		if ((this_M<eps60) .and. (this_M>-eps60)) then ! this_M==0
			return
		else
			this_M = 0.0e0
		end if
	end if

end subroutine fix_MN

subroutine micro_proc
  use vars, only: dudt, dvdt, dwdt, tabs, t, pres, rho, nrestart, qv, gamaz, &
                  qcl, qpl, esatw, qsatw
  use grid, only: nc, dx, dy, dz, dt, rank, nsubdomains, nstep, icycle, z, &
                  dimx1_w, dimx2_w, dimy1_w, dimy2_w
  use params

  implicit none

  integer, parameter :: its=1, ite=nx, jts=1, &
                        jte=ny, kts=1, kte=nzm


! local variable
  integer itimestep
  integer :: i,j,k,kr,ikt,m
  integer :: ii ! loop index [L]
  real, dimension(kts:kte) :: rhocgs, pcgs, zcgs
  real, dimension(naerosol) :: naer1d
  real, dimension(ncloud) :: ncld1d,mcld1d
  real :: Nactivated(nx,ny,nzm), Mactivated(nx,ny,nzm)
  real :: Ncondensation(nx,ny,nzm), Mcondensation(nx,ny,nzm)
  real :: Nanew, Manew, Ncnew, Mcnew
  real :: tt, qq, pp, roro
  real :: dtcond ! time step for condensational growth
  real :: esat, ss_env
  real :: sat_temp
  integer :: bin_begin, bin_end ! speed up growth and collision [L]
  real :: temp_M, temp_N, temp_MN
  real, dimension(ncloud) :: temp_Mall
  real, dimension(ncloud) :: m_ave
  real :: act_temp, N_before_temp, N_after_temp
  real, dimension(ncloud) :: rcld_temp

  itimestep = nstep
 
  if (itimestep.eq.1.and.icycle.eq.1.and.masterproc)then
    if (iceprocs.eq.1)print*,'ICE PROCESES ACTIVE'
    if (iceprocs.eq.0)print*,'LIQUID PROCESES ONLY'
  end if

! set time step
  if (dt<ncond) then
    dtcond = dt
    ncond = 1
  else
    dtcond = 1.0
    ncond = nint(dt/1.0)
  end if

  if (abs(dtcond*ncond-dt)>0.01) then
    print*,"Error: dt is not integer"
    stop
  end if

! change unit to cgs
  do k = kts, kte
    pcgs(k) = pres(k) * 1000.
    rhocgs(k) = rho(k) * 0.001
    zcgs(k) = z(k) * 100.
  end do 



  if (itimestep.eq.1.and.icycle.eq.1) then
    do k = kts, kte
    do j = jts, jte
    do i = its, ite
      do kr = 1, naerosol
        Naer(i,j,k,kr) = numaer0(kr)/rho(k)*1.0e6
      end do
      do kr = 1, ncloud
        Ncld(i,j,k,kr) = numcld0(kr)/rho(k)*1.0e6
        Mcld(i,j,k,kr) = numcld0(kr)*(mcloud(kr)+mcloud2(kr))/2.0/rho(k)*1.0e3
      end do
    end do
    end do
    end do

  else

    do k = kts, kte
    do j = jts, jte
    do i = its, ite
      do kr = 1, ncloud
	    temp_Mall = Mcld(i,j,k,:)
		call fix_MN(Mcld(i,j,k,kr), Ncld(i,j,k,kr), rho(k), kr)
      end do
    end do
    end do
    end do

    do k = kts, kte
    do j = jts, jte
    do i = its, ite

      ! Update Naer
      do kr = 1, naerosol
        if (kr==22) then
          if (inject_aer.and.itimestep.ge.inject_aer_step) then
            Naer(i,j,k,kr) = Naer(i,j,k,kr) + aer_source * dt / rho(k) * 1.0e6
          else
            Naer(i,j,k,kr) = 0.0
          end if
        else
          Naer(i,j,k,kr) = 0.0
        end if
      end do

    end do 
    end do
    end do

  end if ! end of initial step and others

  do k = kts, kte
  do j = jts, jte
  do i = its, ite
    do kr = 1, ncloud
      ncld1d(kr) = Ncld(i,j,k,kr)*rho(k)*1.0e-6
      mcld1d(kr) = Mcld(i,j,k,kr)*rho(k)*1.0e-3
      naer1d(kr) = Naer(i,j,k,kr)*rho(k)*1.0e-6
    end do
! update qc, qr, qna, qnc, qnr
    qh(i,j,k) = sum(mcld1d(1:kcdrop-1))/rho(k)*1.0e3
    qc(i,j,k) = sum(mcld1d(1:krdrop))/rho(k)*1.0e3
    qr(i,j,k) = sum(mcld1d(krdrop+1:ncloud))/rho(k)*1.0e3
	qlwc(i,j,k) = qh(i,j,k) + qc(i,j,k) + qr(i,j,k)
    qna(i,j,k) = sum(naer1d)/rho(k)
	qnh(i,j,k) = sum(ncld1d(1:kcdrop-1))/rho(k)
    qnc(i,j,k) = sum(ncld1d(kcdrop:krdrop))/rho(k)
    qnr(i,j,k) = sum(ncld1d(krdrop+1:ncloud))/rho(k)
	qnt(i,j,k) = qnh(i,j,k) + qnc(i,j,k) + qnr(i,j,k)
    Nactivated(i,j,k) = 0.0
    Mactivated(i,j,k) = 0.0
    Ncondensation(i,j,k) = 0.0
    Mcondensation(i,j,k) = 0.0
    Nanew = 0.0      ! number concentration of newly activated droplets
    Manew = 0.0      ! mass concentration of newly activated droplets
    Ncnew = 0.0      ! number concentration of deactivated aerosol
    Mcnew = 0.0      ! mass concentration change due to condensation

    tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k) + qr(i,j,k)) ! [k]
    qv(i,j,k) = qt(i,j,k) -(qc(i,j,k)+qr(i,j,k)) ! [kg/kg]

    tt = tabs(i,j,k)  ! temperature in [K]
    qq = qv(i,j,k) ! water vapor mixing ratio [kg/kg]
    pp = pcgs(k)       ! pressure
    roro = rhocgs(k) ! density [g/cm3]

    ! Do deliquescence
    !call aerosol_deliq(tt,qq,pp,roro,naer1d,ncld1d,mcld1d,Nanew,Manew)
    call aerosol_deliq2(tt,qq,pp,roro,naer1d,ncld1d,mcld1d,Nanew,Manew)

    tt = t(i,j,k) - gamaz(k) + fac_cond*sum(mcld1d)/rho(k)*1.0e3
    qq = qt(i,j,k) - sum(mcld1d)/rho(k)*1.0e3

    ! Judge growth and collision bin range [L]
    bin_begin=1
    bin_end=ncloud
    do ii=1,ncloud
        if (ncld1d(ii).gt.eps60) then
            bin_begin=ii
            EXIT
        end if
    end do
    do ii=ncloud,bin_begin,-1
        if (ncld1d(ii).gt.eps60) then
            bin_end=ii
            EXIT
        else if (ii==bin_begin) then
            bin_end = bin_begin
        end if
    end do

    N_before_temp = sum(ncld1d(kcdrop:ncloud))

    ! Do growth
    do ikt = 1, ncond
      call cloud_cond_growth(tt,qq,pp,roro,dtcond,ncld1d,mcld1d,Ncnew,bin_begin, bin_end,i,j,k)
      if (bin_begin .gt. bin_end) then
        print*, 'WARNING: bin_begin>bin_end after cloud_cond_growth',bin_begin, bin_end
        bin_end = bin_begin
      end if
      Ncondensation(i,j,k) = Ncondensation(i,j,k) + Ncnew ! Ncnew: evaporated N
    end do

    ! statistics for activation and deactivation
    N_after_temp = sum(ncld1d(kcdrop:ncloud))
	act_temp = N_after_temp - N_before_temp
    if (act_temp>0e0) then
	  act_rate(i,j,k) = act_temp/rho(k)/dt
	  deact_rate(i,j,k) = 0e0
	else
	  act_rate(i,j,k) = 0e0
	  deact_rate(i,j,k) = -act_temp/rho(k)/dt
	end if

    ! Do collision
    if (docoll) call cloud_collision(dt, ncld1d, mcld1d,bin_begin, bin_end) ! [L]

    ! Update results to 3D array
    do kr = 1, naerosol
       if (kr==22) then
         Naer(i,j,k,kr) = naer1d(kr)/rho(k)*1e6
       else
         Naer(i,j,k,kr) = naer1d(kr)/rho(k)*1e6
       end if
    end do

    do kr = 1, ncloud
      Ncld(i,j,k,kr) = ncld1d(kr)/rho(k)*1e6
      Mcld(i,j,k,kr) = mcld1d(kr)/rho(k)*1e3
	  if (Ncld(i,j,k,kr) .le. 1e-10) then
	    rcld_temp(kr) = 0e0
	  else
	    rcld_temp(kr) = (3.0/4.0/pi*mcld1d(kr)/ncld1d(kr))**(1.0/3.0)
	  end if
    end do

    qh(i,j,k) = sum(mcld1d(1:kcdrop-1))/rho(k)*1.0e3
    qc(i,j,k) = sum(mcld1d(1:krdrop))/rho(k)*1.0e3
    qr(i,j,k) = sum(mcld1d(krdrop+1:ncloud))/rho(k)*1.0e3
	qlwc(i,j,k) = qh(i,j,k) + qc(i,j,k) + qr(i,j,k)
    qna(i,j,k) = sum(naer1d)/rho(k)
	qnh(i,j,k) = sum(ncld1d(1:kcdrop-1))/rho(k)
    qnc(i,j,k) = sum(ncld1d(kcdrop:krdrop))/rho(k)
    qnr(i,j,k) = sum(ncld1d(krdrop+1:ncloud))/rho(k)
	qnt(i,j,k) = qnh(i,j,k) + qnc(i,j,k) + qnr(i,j,k)

    if (qv(i,j,k) .le. 1.e-15) qv(i,j,k) = 1.e-15
	if (qh(i,j,k) .le. 1.e-15) qh(i,j,k) = 0e0
    if (qc(i,j,k) .le. 1.e-15) qc(i,j,k) = 0e0
    if (qr(i,j,k) .le. 1.e-15) qr(i,j,k) = 0e0
	if (qlwc(i,j,k) .le. 1.e-15) qlwc(i,j,k) = 0e0

    tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond*(qc(i,j,k)+qr(i,j,k))

! micro diagnose
    qv(i,j,k) = qt(i,j,k) - (qc(i,j,k)+qr(i,j,k))
    qcl(i,j,k) = qc(i,j,k) + qr(i,j,k)
    qpl(i,j,k) = 0.0

    sat_temp = qsatw(tabs(i,j,k),pp*1.e-3)
    sat_temp = qv(i,j,k)/sat_temp - 1.0
    ssatw(i,j,k) = sat_temp*1.0e2


    if (sum(ncld1d)>1.0e-15) then
       rcmean(i,j,k) = sum(rcld_temp(kcdrop:ncloud)*ncld1d(kcdrop:ncloud))/ &
	                   sum(ncld1d(kcdrop:ncloud))
       rcsig(i,j,k) = (sum(ncld1d(kcdrop:ncloud)*(rcld_temp(kcdrop:ncloud) - & 
	                  rcmean(i,j,k))**2)/sum(ncld1d(kcdrop:ncloud)))**0.5
       rceff(i,j,k) = sum(rcld_temp(kcdrop:ncloud)**3*ncld1d(kcdrop:ncloud)) / &
	                  sum(rcld_temp(kcdrop:ncloud)**2*ncld1d(kcdrop:ncloud))
       rel_dis(i,j,k) = rcsig(i,j,k)/rcmean(i,j,k)
    else
       rcmean(i,j,k) = 0.0
       rcsig(i,j,k) = 0.0
       rceff(i,j,k) = 0.0
       rel_dis(i,j,k) = 0.0
    end if
  end do
  end do
  end do

end subroutine micro_proc

subroutine micro_statistics
  use vars
  use hbuffer, only: hbuf_put, hbuf_avg_put
  use params

  real :: tmp(2), factor_xy

! bulk properties
  real, dimension(nzm) :: na_b, nh_b, nc_b, nt_b
  real, dimension(nzm) :: qh_b, qc_b, qt_b
  real, dimension(nzm) :: rmean_b, rsig_b, reff_b, reldis_b
  real, dimension(nzm) :: act_b, deact_b
  integer :: nx_b, ny_b, offset_f, offset_b
  integer :: xsb, xeb, ysb, yeb, zsb, zeb
  real :: factor_xy_b
  real :: factor_vol(nzm)
  integer :: i,j,k,m

  factor_xy = 1.0 / float(nx*ny)

  offset_f = 0
  offset_b = 0
  na_b = 0e0
  nh_b = 0e0
  nc_b = 0e0
  nt_b = 0e0
  qh_b = 0e0
  qc_b = 0e0
  qt_b = 0e0
  rmean_b = 0e0
  rsig_b = 0e0
  reff_b = 0e0
  reldis_b = 0e0
  act_b = 0e0
  deact_b = 0e0

  xsb = 1
  xeb = nx
  ysb = 1
  yeb = ny
  zsb = 1 + offset_f
  zeb = nzm - offset_b

    if(mod(rank,nsubdomains_x).eq.0) then 
        xsb = 1 + offset_f
    end if
    if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
        xeb = max(xsb, xeb - offset_b)
    end if   
    if(rank.lt.nsubdomains_x) then
        ysb = 1 + offset_f
    end if   
    if(rank.gt.nsubdomains-nsubdomains_x-1) then
        yeb = max(ysb,yeb - offset_b)
    end if   
    if(xsb.eq.xeb) then
        print*,"domain too small, reset offset_f"
        stop
        zsb=zeb+2
        ysb=yeb+2
        xsb=xeb+2
    end if
    if(ysb.eq.yeb) then
        print*, "domain too small, reset offset_f"
        stop
        zsb=zeb+2
        ysb=yeb+2
        xsb=xeb+2
    end if
    nx_b = nx_gl - offset_b - offset_f
    ny_b = ny_gl - offset_b - offset_f
    factor_xy_b = float(nsubdomains_x*nsubdomains_y)/float(nx_b*ny_b)

    do k=zsb,zeb
      do i=xsb,xeb
        do j=ysb,yeb
          na_b(k) = na_b(k) + qna(i,j,k)
		  nh_b(k) = nh_b(k) + qnh(i,j,k)
          nc_b(k) = nc_b(k) + qnc(i,j,k) + qnr(i,j,k)
		  nt_b(k) = nt_b(k) + qnt(i,j,k)
		  qh_b(k) = qh_b(k) + qh(i,j,k)
		  qc_b(k) = qc_b(k) + qc(i,j,k) + qr(i,j,k)
		  qt_b(k) = qt_b(k) + qlwc(i,j,k)
          rmean_b(k) = rmean_b(k) + rcmean(i,j,k)
          rsig_b(k) = rsig_b(k) + rcsig(i,j,k)
          reff_b(k) = reff_b(k) + rceff(i,j,k)
          reldis_b(k) = reldis_b(k) + rel_dis(i,j,k)
		  act_b(k) = act_b(k) + act_rate(i,j,k)
		  deact_b(k) = deact_b(k) + deact_rate(i,j,k)
        end do
      end do
    end do

    call hbuf_put('BNA',na_b,factor_xy_b)
    call hbuf_put('BNH',nh_b,factor_xy_b)
    call hbuf_put('BNC',nc_b,factor_xy_b)
    call hbuf_put('BNT',nt_b,factor_xy_b)
    call hbuf_put('BQCh',qh_b,factor_xy_b)
    call hbuf_put('BQCl',qc_b,factor_xy_b)
    call hbuf_put('BQCt',qt_b,factor_xy_b)
    call hbuf_put('BRMEAN',rmean_b,1.0e4*factor_xy_b)
    call hbuf_put('BRSIG',rsig_b,1.0e4*factor_xy_b)
    call hbuf_put('BREFF',reff_b,1.0e4*factor_xy_b)
    call hbuf_put('BRELDIS',reldis_b,factor_xy_b)
	call hbuf_put('BARATE',act_b,factor_xy_b)
	call hbuf_put('BDRATE',deact_b,factor_xy_b)

end subroutine micro_statistics

subroutine micro_print
end subroutine micro_print


subroutine micro_hbuf_init(namelist, deflist, unitlist, status, average_type, count, trcount)
  character(*) namelist(*), deflist(*), unitlist(*)
  integer status(*), average_type(*), count, trcount

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QTFLUX'
  deflist(count) = 'Nonprecipitating water flux (Total)'
  unitlist(count) = 'W/m2'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QTFLUXS'
  deflist(count) = 'Nonprecipitating water flux (SGS)'
  unitlist(count) = 'W/m2'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QPFLUX'
  deflist(count) = 'Precipitating-water turbulent flux (Total)'
  unitlist(count) = 'W/m2'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QPFLUXS'
  deflist(count) = 'Precipitating-water turbulent flux (SGS)'
  unitlist(count) = 'W/m2'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QC'
  deflist(count) = 'Cloud water (microphysics)'
  unitlist(count) = 'g/m^3'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QR'
  deflist(count) = 'Rain water (microphysics)'
  unitlist(count) = 'g/m^3'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'SSW'
  deflist(count) = 'Supersaturation with respect to water'
  unitlist(count) = '%'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'NA'
  deflist(count) = 'Aerosol number concentration'
  unitlist(count) = 'cm-3'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'NC'
  deflist(count) = 'Cloud droplet number concentration'
  unitlist(count) = 'cm-3'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'NR'
  deflist(count) = 'Rain drop number concentration'
  unitlist(count) = 'cm-3'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QLCODTEND'
  deflist(count) = 'Liquid mixing ratio tendency due to cond/evap'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QLCOLTEND'
  deflist(count) = 'Liquid mixing ratio tendency due to collision'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QLSEDLTEND'
  deflist(count) = 'Liquid mixing ratio tendency due to sed'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QVADVT'
  deflist(count) = 'Water vapor tendency due to advection'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QLADVT'
  deflist(count) = 'Liquid water tendency due to advection'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QVTENDT'
  deflist(count) = 'Total water vapor tendency'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'QLTENDT'
  deflist(count) = 'Total liquid water tendency'
  unitlist(count) = 'g/kg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'REFFC'
  deflist(count) = 'Effective Radius'
  unitlist(count) = 'um'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'RMEAN'
  deflist(count) = 'Mean Radius'
  unitlist(count) = 'um'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'NACT'
  deflist(count) = 'Number of activated droplets'
  unitlist(count) = '#/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'NLOSS'
  deflist(count) = 'Droplet loss due to collision coalescence'
  unitlist(count) = '#/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'NDSEDTEND'
  deflist(count) = 'Droplet number tendency due to sed'
  unitlist(count) = '/cm3/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'DLTot'
  deflist(count) = 'Total diffusion growth'
  unitlist(count) = 'g/kg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'CLTot'
  deflist(count) = 'Total collisions'
  unitlist(count) = '#/cm3'
  status(count) = 1
  average_type(count) = 0

! add for Bulk properties
! Fan Yang 11/21/2020

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BNA'
  deflist(count) = 'Bulk aerosol number concentration'
  unitlist(count) = '#/mg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BNH'
  deflist(count) = 'Bulk haze number concentration'
  unitlist(count) = '#/mg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BNC'
  deflist(count) = 'Bulk cloud droplet number concentration'
  unitlist(count) = '#/mg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BNT'
  deflist(count) = 'Bulk total particle number concentration'
  unitlist(count) = '#/mg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BQCh'
  deflist(count) = 'Bulk haze content'
  unitlist(count) = 'g/kg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BQCl'
  deflist(count) = 'Bulk cloud content'
  unitlist(count) = 'g/kg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BQCt'
  deflist(count) = 'Bulk total water content'
  unitlist(count) = 'g/kg'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BRMEAN'
  deflist(count) = 'Bulk mean cloud droplet radius'
  unitlist(count) = 'um'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BRSIG'
  deflist(count) = 'Bulk standard deviation of cloud droplet size distribution'
  unitlist(count) = 'um'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BREFF'
  deflist(count) = 'Bulk effective radius'
  unitlist(count) = 'um'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BRELDIS'
  deflist(count) = 'Bulk relative dispersion'
  unitlist(count) = ' '
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BARATE'
  deflist(count) = 'Bulk activation rate'
  unitlist(count) = '#/mg/s'
  status(count) = 1
  average_type(count) = 0

  count = count + 1
  trcount = trcount + 1
  namelist(count) = 'BDRATE'
  deflist(count) = 'Bulk deactivation rate'
  unitlist(count) = '#/mg/s'
  status(count) = 1
  average_type(count) = 0

endsubroutine micro_hbuf_init

subroutine micro_flux
use vars, only: fluxbq, fluxtq, fluxlq, fluxrq, fluxqq, fluxhq

  fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
  fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
  fluxlmk(:,:,:) = 0. ! initialize all fluxes at left wall to zero
  fluxrmk(:,:,:) = 0. ! initialize all fluxes at right wall of domain to zero
  fluxqmk(:,:,:) = 0. ! initialize all fluxes at front wall to zero
  fluxhmk(:,:,:) = 0. ! initialize all fluxes at back wall of domain to zero
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)
  fluxlmk(:,:,index_water_vapor) = fluxlq(:,:)
  fluxrmk(:,:,index_water_vapor) = fluxrq(:,:)
  fluxqmk(:,:,index_water_vapor) = fluxqq(:,:)
  fluxhmk(:,:,index_water_vapor) = fluxhq(:,:)

end subroutine micro_flux

subroutine micro_diagnose
  use vars
  integer i,j,k

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        qv(i,j,k) = qt(i,j,k) - (qc(i,j,k) + qr(i,j,k))
        qcl(i,j,k) = qc(i,j,k) + qr(i,j,k)
        qci(i,j,k) = 0.0
        qpl(i,j,k) = 0.0
        qpi(i,j,k) = 0.0
      end do
    end do
  end do
end subroutine micro_diagnose

!-----------------------------------------------------------------------
! Function that computes total water in a domain:
! Don't change this one.

real function total_water() 

  use vars, only : nstep,nprint,adz,dz,rho

  integer k,m

  total_water = 0. 
  if(mod(nstep,nprint).ne.0) return

  do m=1,nmicro_fields

   if(flag_wmass(m).eq.1) then 

    do k=1,nzm
      total_water = total_water + &
       sum(micro_field(1:nx,1:ny,k,m))*adz(k)*dz*rho(k)
    end do

   end if

  end do

end function total_water

subroutine micro_precip_fall()
  use vars
  integer :: hydro_type
  real :: omega(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  integer :: i, j, k, m, kr
  integer, parameter :: its=1, ite=nx, jts=1,     &    
                         jte=ny, kts=1, kte=nzm    

  real,  DIMENSION(kts:kte):: rhocgs

! Initialize arrays that accumulate surface precipitation flux
                                                                                         
 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then 
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if


 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
 end do

! update
  do k = kts,kte
  do j = jts,jte
  do i = its,ite
    qc(i,j,k) = 0.0
    qr(i,j,k) = 0.0
    do m = 1, ncloud
       if (m.le.krdrop) then
         qc(i,j,k) = qc(i,j,k) + Mcld(i,j,k,m)
       else
         qr(i,j,k) = qr(i,j,k) + Mcld(i,j,k,m)
       end if
    end do
    if (qc(i,j,k).lt.1.e-15) qc(i,j,k) = 0.0
    if (qr(i,j,k).lt.1.e-15) qr(i,j,k) = 0.0
    qv(i,j,k) = qt(i,j,k) - (qc(i,j,k) + qr(i,j,k))
    tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k) + qr(i,j,k))
  end do
  end do
  end do

  hydro_type = 0 ! for mass

  do m = 1,ncloud
    omega(:,:,:) = Mcld(:,:,:,m) 
    call precip_fall(omega, term_vel_cloud, hydro_type, omega, m)
    Mcld(:,:,:,m) = omega(:,:,:)
  end do

  hydro_type = 3 ! for number

  do m = 1,ncloud
    omega(:,:,:) = Ncld(:,:,:,m)
	call precip_fall(omega, term_vel_cloud, hydro_type, omega, m)
	Ncld(:,:,:,m) = omega(:,:,:)
  end do

! update
  do k = kts,kte
  do j = jts,jte
  do i = its,ite
    qc(i,j,k) = 0.0
    qr(i,j,k) = 0.0
    do m = 1, ncloud
       if (m.le.krdrop) then
         qc(i,j,k) = qc(i,j,k) + Mcld(i,j,k,m)
       else
         qr(i,j,k) = qr(i,j,k) + Mcld(i,j,k,m)
       end if
    end do
    if (qc(i,j,k).lt.1.e-15) qc(i,j,k) = 0.0
    if (qr(i,j,k).lt.1.e-15) qr(i,j,k) = 0.0
    qt(i,j,k) = qv(i,j,k) + qc(i,j,k) + qr(i,j,k)
    tabs(i,j,k) = t(i,j,k) - gamaz(k) + fac_cond * (qc(i,j,k) + qr(i,j,k))
  end do
  end do
  end do

! surface precipitation area fraction statistics
  do j=1,ny
    do i=1,nx
      if((qr(i,j,1)+qc(i,j,1)).gt.1.e-6) s_ar = s_ar+dtfactor
    end do
  end do
end subroutine micro_precip_fall

real function term_vel_cloud(i,j,k,nbin)

  real, intent(in) :: i,j,k
  integer :: nbin

  term_vel_cloud = cloud_term_vel(nbin) *0.01

end function term_vel_cloud

! v6.9.4
function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi

function Get_rmean() ! liquid water mean
  real, dimension(nx,ny,nzm) :: Get_rmean
  Get_rmean = rmean
end function Get_rmean


end module
