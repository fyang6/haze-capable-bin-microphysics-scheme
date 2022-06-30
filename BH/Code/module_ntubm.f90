module module_ntubm

! This module is a double-moment bin microphysical schemes based on Chen and
! Lamb JAS 1994
! Reference:
! Chen, J.P. and Lamb, D., 1994. Simulation of cloud microphysical and chemical
! processes using a multicomponent framework. Part I: Description of the
! microphysical model. Journal of the Atmospheric Sciences, 51(18),
! pp.2613-2630.
!
! The code is implemented in SAM by Fan Yang @ BNL (fanyang@bnl.gov), 2021
!
! History:
! version, date, description
! v0.1(F. Yang); 01/12/2021; only warm cloud microphysical processes are included
! v0.2(F. Yang); 06/21/2021; bug fixed
! v0.3(F. Yang); 07/03/2021; this version includes haze particles, growth equation updated

use params, only: cp, ggr, rgas, rv, fac_cond, lcond, pi, diffelq, therco 
use grid, only: masterproc, rank
use vars, only: esatw, qsatw
use micro_prm

implicit none

real, dimension(naerosol), save :: raerosol ! [cm] aerosol radius
real, dimension(naerosol), save :: numaer0 ! [cm-3] initial aerosol number concentration
real, dimension(ncloud), save :: numcld0 ! [cm-3] initial cloud number concentration
real, dimension(ncloud), save :: rcloud ! [cm] cloud droplet radius
real, dimension(ncloud), save :: mcloud ! [g] cloud droplet mass
real, dimension(ncloud), save :: mcloud2 ! edge

real, dimension(ncloud), save :: cloud_term_vel ! [cm/s] terminal velocity of cloud droplets

real, dimension(ncloud,ncloud), save :: kernel_cloud ! input array of kernels for cloud droplet collision
logical, save :: bug_happen

contains
        
subroutine cloud_collision(dtcond, ncld, mcld, bin_begin, bin_end) ! [L]
    real, INTENT(IN) :: dtcond
    integer, INTENT(IN) :: bin_begin, bin_end
    real, dimension(ncloud), INTENT(INOUT) :: ncld, mcld
    real, dimension(ncloud) :: E1A, E2A, deltaM, deltaN
    integer :: i, ii, jj
    real :: m1, m2, Ek, Em, En, dMij, dNij, dMi, dMj
    E1A = mcloud
    E2A = mcloud2
    deltaM = 0e0
    deltaN = 0e0

    do i = bin_begin, bin_end
        call get_linear_eq_E(E1A(i), E2A(i), mcld(i), ncld(i), Ek, Em, En)
        if (bug_happen) then
            print *, 'bug_happen: in cloud_collision begin, bin',i
            call sleep(2)
            stop
        end if
    end do

    do ii = bin_begin, bin_end
       do jj = ii, bin_end
            if (ncld(ii)>epsS .and. ncld(jj)>epsS) then
                call delta_MN(ii,jj,mcld(ii), ncld(ii), mcld(jj), ncld(jj), dtcond, dMij, dNij, dMi, dMj)
                
                if (dNij > 0e0) then
                    m1 = E1A(ii) + E1A(jj)
                    m2 = E2A(ii) + E2A(jj)
                    deltaM(ii) = deltaM(ii) - dMi
                    deltaN(ii) = deltaN(ii) - dNij
                    deltaM(jj) = deltaM(jj) - dMj
                    deltaN(jj) = deltaN(jj) - dNij

                    if (m2 .lt. mcloud2(jj)) then
                        deltaM(jj) = deltaM(jj) + dMij
                        deltaN(jj) = deltaN(jj) + dNij
                    else
                        call get_linear_eq_E(m1, m2, dMij, dNij, Ek, Em, En)
                        if (bug_happen) then
                            print *, 'bug_happen: error in cloud_collision collided bin',ii,'x',jj
                            call sleep(2)
                            stop
                        end if
                        if (m2 .lt. mcloud2(jj)) then
                            deltaM(jj) = deltaM(jj) + dMij
                            deltaN(jj) = deltaN(jj) + dNij
                        else
                            call redistribute_cloud(m1, m2, dMij, dNij, Ek, Em, En, deltaM, deltaN, jj)
                        end if
                    end if
                end if
            end if
        end do
    end do

    mcld = mcld + deltaM
    ncld = ncld + deltaN
    
end subroutine cloud_collision

subroutine delta_MN(ii,jj,M1,N1,M2,N2,delta_t,dM, dN, dMi, dMj)
     integer, INTENT(IN) :: ii, jj
     real, INTENT(IN) :: M1,N1,M2,N2,delta_t
     real, INTENT(OUT) :: dM, dN, dMi, dMj
     real :: K, r1, r2, K2

     K = kernel_cloud(ii,jj)

     dN = N1 * N2 * K * delta_t
     if (ii==jj) then
         dN = dN/2e0
     end if
     dMi = dN * (M1/N1)
     dMj = dN * (M2/N2)
     dM = dN * (M1/N1 + M2/N2)
end subroutine delta_MN

subroutine r2K(r1,r2,K)
    real, INTENT(IN) :: r1, r2
    real, INTENT(OUT) :: K
    real :: vl, vr, E

    call long_kernel_efficiency(r1, r2, E)
    
    call calc_terminal_velocity(r1, vl)
    call calc_terminal_velocity(r2, vr)
    
    K = pi * (r1 + r2)**2e0 * E * abs(vr-vl)
end subroutine r2K

subroutine long_kernel_efficiency(r1, r2, E)
     real, INTENT(IN) :: r1, r2
     real, INTENT(OUT) :: E
     real :: thres, rL, rr
     thres = 50e-4
     
     call r1_r2(r1,r2,rL,rr)

     if (rr >= thres) then
         E = 1e0
     else if (rr <= 0 ) then

         print *, "Error: r<0 in long_kernel_efficiency "
         call exit(0)
     else
         E = 4.5e4 * rr**2e0 * (1e0 -3e-4/rL)
         if (E<10e-3) then
                 E = 10e-3
         end if
     end if
end subroutine long_kernel_efficiency

subroutine r1_r2(r1,r2,rL,rr)
   real, INTENT(IN) :: r1, r2
   real, INTENT(OUT) :: rL, rr
     if (r1 > r2) then
             rr = r1
             rL = r2
     else
             rr = r2
             rL = r1
     end if
end subroutine r1_r2

subroutine calc_terminal_velocity(r,v)
     real, INTENT(IN) :: r
     real, INTENT(OUT) :: v
     
     if (r>0.7e0) then
         call terminal_velocity_1070_7000(r,v)
     else if (r >=0.107e0) then
         call terminal_velocity_1070_7000(r,v)
     else if (r >=0.0019e0) then
         call terminal_velocity_19_1070(r,v)
     else if (r >= 0.5e-4) then
         call terminal_velocity_05_19(r,v)
     else
         call terminal_velocity_05_19(r,v)
     end if
end subroutine calc_terminal_velocity

subroutine polynomial(X,b_list,Y)
      real, intent(in) :: X, b_list(:)
      real, intent(out) :: Y
      integer :: n, i

      n = size(b_list)
      Y = 0e0
      do i = 1, n
         Y = Y + X**(i-1) * b_list(i)
      end do
end subroutine polynomial

subroutine terminal_velocity_1070_7000(r,v)
     real, INTENT(IN) :: r
     real, INTENT(OUT) :: v
     real :: b_list(6), C3, B0, Np, X, Y, Nre
	 real :: sigma_H2O
	 sigma_H2O = 72.8
     b_list = [-0.500015e1, 0.523778e1, -0.204914e1, 0.475294e0, -0.542819e-1, 0.238449e-2]
        
     C3 = 4e0 * delta_rho * G0 /3e0 / sigma_H2O
     B0 = C3 * r**2
     Np = sigma_H2O**3 * rho_air**2 / eta0**4 / delta_rho / G0
     X = log(B0 * Np**(1e0/6e0))
     
     call polynomial(X, b_list, Y)
     Nre = Np**(1e0/6e0) * exp(Y)

     v = eta0 * Nre / rho_air / r
end subroutine terminal_velocity_1070_7000

subroutine calc_Csc(r, Csc)
     real, INTENT(IN) :: r
     real, INTENT(OUT) :: Csc
     Csc = 1e0 + 2.51e0 * L0 / r    
end subroutine calc_Csc

subroutine terminal_velocity_19_1070(r,v)
     real, INTENT(IN) :: r
     real, INTENT(OUT) :: v
     real :: b_list(7), C2, Nda, X, Y, Csc, Nre
     b_list = [-0.318657e1, 0.992696e0, -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4, -0.327815e-5]
        
     C2 = 4e0 * rho_air * delta_rho * G0 /3e0 / eta0**2
     Nda = C2 * r**3
     X = log(Nda)
     
     call polynomial(X, b_list, Y)
     call calc_Csc(r, Csc)
     Nre = Csc * exp(Y)

     v = eta0 * Nre / rho_air / r
end subroutine terminal_velocity_19_1070


subroutine terminal_velocity_05_19(r,v)
     real, INTENT(IN) :: r
     real, INTENT(OUT) :: v
     real :: C1, Csc
        
     call calc_Csc(r, Csc)
     C1 = delta_rho * G0 / 18e0 / eta0

     v = C1 * Csc * r**2
end subroutine terminal_velocity_05_19

subroutine cloud_cond_growth(tt, qq, pp, roro, dtcond, ncld, mcld, ncnew, bin_begin, bin_end,ii,jj,kk)
    real, INTENT(IN) :: pp, roro, dtcond
	integer, INTENT(IN), optional :: ii,jj,kk
    integer, INTENT(INOUT) :: bin_begin, bin_end
    real, INTENT(INOUT) :: tt, qq
    real, INTENT(OUT) :: ncnew
    real, dimension(ncloud), INTENT(INOUT) :: ncld, mcld
    real, dimension(ncloud) :: Nnew, Mnew
    real :: Gcond, ss_env
	real :: A_factor, B_factor
    real, dimension(ncloud) :: N_temp
    real, dimension(ncloud) :: M_N_temp
    real :: M_Ni, m1_temp,m2_temp, Em, Ek, En
    integer :: i

    ncnew = 0.0

    Nnew = 0.0
    Mnew = 0.0

    call growth_factor(tt,qq,pp,Gcond,ss_env,A_factor,B_factor,ii,jj,kk)
    do i = bin_begin, bin_end
        if (ncld(i) > 0 ) then

            call do_one_cloud_growth(i, mcld(i), ncld(i), Mnew, Nnew, ncnew, Gcond, &
			                         A_factor, B_factor, ss_env, dtcond,ii,jj,kk)
        end if 
    end do
    call get_new_bin_begin(Nnew,bin_begin)
    call get_new_bin_end(Nnew,bin_end)

    ncld = Nnew
    mcld = Mnew
end subroutine cloud_cond_growth

subroutine do_one_cloud_growth(i, M1, N1, M2, N2, N3, G, Af, Bf, s_1, dt,ii,jj,kk)
	integer, INTENT(IN), optional :: ii,jj,kk
    integer, INTENT(IN) :: i
    real, INTENT(IN) :: G, s_1, dt, Af, Bf
    real, INTENT(INOUT) :: M1, N1
    real, dimension(ncloud), INTENT(INOUT):: M2, N2
    real, INTENT(INOUT) :: N3
    real :: mm1, mm2, Mn, Ek, Em, En
    real :: ss
    logical :: debug = .False.
	integer :: this_ind
	real :: rrr1,rrr2,rrr3,rrr4
	real :: mm1_old, mm2_old

    ss = s_1 - 1.0

    mm1 = mcloud(i)
    mm2 = mcloud2(i)
    
    if (M1<eps60) then
        M1 = 0e0
        N1 = 0e0
        RETURN
    else if (N1<eps) then 
        M2(i) = M2(i) + M1
        N2(i) = N2(i) + N1
        RETURN
    end if

    call get_linear_eq_E(mm1, mm2, M1, N1, Ek, Em, En)

    if (bug_happen) then
        print*, 'bug_happen: after 1st linear in do_one'
        call sleep(2)
        stop
    end if
    
	mm1_old = mm1
	mm2_old = mm2
    call m_end_grow(mm1, G, ss, Af, Bf, dt,ii,jj,kk) 
    call m_end_grow(mm2, G, ss, Af, Bf, dt,ii,jj,kk) 
	if (mm1<=mcloud(1)) then
		mm1 = mcloud(1) + epsS
	end if
	if (mm2<= mm1) then
		mm2 = mm1 + epsS 
	end if
   call M_sum_grow(M1, Mn, N1, G, ss, Af, Bf, dt, mm1, mm2,ii,jj,kk)

   if (abs(mm2-mm1)<epsL) then
		call locate_m((mm2+mm1)/2.0e0,this_ind) 
		M2(this_ind) = M2(this_ind) + Mn
		N2(this_ind) = N2(this_ind) + N1
	else
	   call get_linear_eq_E(mm1, mm2, Mn, N1, Ek, Em, En) 
	   if (bug_happen) then
		  print*, 'bug_happen: after 2nd linear in do_one'
		  call sleep(2)
		  stop
	   end if
	   call redistribute_cloud(mm1, mm2, Mn, N1, Ek, Em, En, M2, N2, i)
   end if

end subroutine do_one_cloud_growth

    
subroutine locate_m(mm,ind1,if_warning)
    real, INTENT(IN) ::mm
    integer, INTENT(OUT) :: ind1
    logical, optional :: if_warning
    logical :: print_warning
    integer :: i
    
    if (present(if_warning)) then
        print_warning = if_warning
    else
        print_warning = .TRUE.
    end if

    do i = 1, ncloud
        if (mcloud(i)<=mm .and. mcloud2(i)>=mm) then
            ind1 = i
            
            RETURN
        end if
    end do

    if (print_warning) then
        print*, 'WARNING: locate_m not found ind for the m',mm
    end if
    ind1 = -999
end subroutine locate_m

subroutine locate_m_new(mm,ind1,left_i,right_i,if_warning)
    real, INTENT(IN) ::mm
    integer, INTENT(IN) :: left_i, right_i
    integer, INTENT(OUT) :: ind1
    logical, optional :: if_warning
    logical :: print_warning
    integer :: i
    
    if (present(if_warning)) then
        print_warning = if_warning
    else
        print_warning = .TRUE.
    end if

    do i = left_i, right_i
        if (mcloud(i)<=mm .and. mcloud2(i)>=mm) then
            ind1 = i
            RETURN
        end if
    end do

    if (print_warning) then
        print*, 'WARNING: locate_m not found ind for the m',mm
    end if
    ind1 = -999
end subroutine locate_m_new

subroutine redistribute_cloud(mm1_ori, mm2_ori, M, N, Ek, Em, En, M2, N2,guess)
    real, INTENT(IN):: mm1_ori, mm2_ori, M, N, Ek, Em, En
    integer, INTENT(IN) :: guess
    real, dimension(ncloud), INTENT(INOUT) :: M2, N2
    integer :: i ,j,cc, j1, j2
    real ::  div_temp
    real :: m1n, m2n, dM_temp, dN_temp
    real, dimension(10) :: save_M, save_N, save_m1, save_m2, save_bs
    integer,dimension(10) :: save_i
	logical, dimension(10) :: test1, test2, test3, test4, testA
    integer :: left_i, right_i, c_temp
	real :: M_N
	integer :: min_bs_ind_temp, min_bs_ind
	real, dimension(4) :: save_M_temp
	logical :: strict

	M_N = M/N

	if ((M_N>mm2_ori) .or. (M_N<mm1_ori)) then
		print*, 'Error: Test0.1 shift at the beginning of redistribution'
		go to 9999
	end if

    left_i = guess
    right_i = guess
    do while (mm1_ori < mcloud(left_i))
        left_i = left_i - 1
    end do
    do while (mm1_ori > mcloud2(left_i))
        left_i = left_i + 1
    end do
    do while (mm2_ori > mcloud2(right_i))
        right_i = right_i + 1
    end do
    do while (mm2_ori < mcloud(right_i))
        right_i = right_i - 1
    end do

    ! Case 1: in one bin
    if (left_i == right_i) then
        M2(left_i) = M2(left_i) + M
        N2(left_i) = N2(left_i) + N
        RETURN
    else if (left_i > right_i) then
        print*, 'Error: Test0.2 redistribution find wrong ind edge', left_i, right_i, mm1_ori, mm2_ori
		go to 9999
    end if

    ! Case 2: exceed edge, then ignore all 
    if (left_i<1) then
        print*, 'Infor: evaportae', mm1_ori, mm2_ori, M, N
        RETURN
    else if (right_i>ncloud) then
        print *, 'Infor: too large',mm1_ori, mm2_ori, M, N
        RETURN
    end if

    ! Case 3: extreme small bin, or very few N
    if (mm2_ori-mm1_ori<epsL .or. N<eps) then
        call locate_m_new(M_N, j2, left_i, right_i,.true.)
        M2(j2) = M2(j2) + M
        N2(j2) = N2(j2) + N
        RETURN
    end if 

    ! Case 4: distribute to left_i ~ right_i bins
	save_m1 = 0.0e0
	save_m2 = 0.0e0
	save_bs = 0.0e0
	save_i = 0
	save_M = 0.0e0
	save_N = 0.0e0

	test3 = .FALSE.
	test4 = .FALSE.

    cc = 1 ! use 1 as an empty parameter
    do i = left_i, right_i
        cc = cc + 1
        if (i == left_i) then
            save_m1(cc) = mm1_ori
        else
            save_m1(cc) = mcloud(i)
        end if
        if (i == right_i) then
            save_m2(cc) = mm2_ori
        else
            save_m2(cc) = mcloud2(i)
        end if
        save_bs(cc) = save_m2(cc) - save_m1(cc)
        save_i(cc) = i
    end do

	if (Ek>0) then
		min_bs_ind = 2
	else
		min_bs_ind = cc
	end if
	
	do i = 2, cc
		if (i .ne. min_bs_ind) then
			call line2MN(save_m1(i), save_m2(i), Ek, Em, En, save_M(i), save_N(i))
		end if
	end do
	save_M(min_bs_ind) = M - sum(save_M(2:cc))
	save_N(min_bs_ind) = N - sum(save_N(2:cc))

	
	min_bs_ind_temp = minval(minloc(save_N(2:cc))) + 1
	if (min_bs_ind_temp .ne. min_bs_ind) then
		call line2MN(save_m1(min_bs_ind), save_m2(min_bs_ind), Ek, Em, En, &
					save_M(min_bs_ind), save_N(min_bs_ind))
		min_bs_ind = min_bs_ind_temp
		save_M(min_bs_ind) = 0.0e0
		save_N(min_bs_ind) = 0.0e0
		save_M(min_bs_ind) = M - sum(save_M(2:cc))
		save_N(min_bs_ind) = N - sum(save_N(2:cc))
	end if

	do j = 2, cc 
		! Test 3: negative value
		if (save_M(j)<0 .or. save_N(j)<0) then

			if (j==2) then
				j1 = 2
				j2 = 3
			else if (j==cc) then
				j1 = cc
				j2 = cc-1
			else
				print*, 'Warning: test3 fix failed',j-1,'in',cc-1
			end if
			save_M(j2) = save_M(j2) + save_M(j1)
			save_N(j2) = save_N(j2) + save_N(j1)
			save_M(j1) = 0.0e0
			save_N(j1) = 0.0e0
		end if 
	end do 

	do j = 2, cc
		! Test 4: shift to other bins
        call locate_m_new(save_M(j)/save_N(j), j1, save_i(j), save_i(j),.FALSE.)
        if (j1 == -999) then 
			test4(j) = .TRUE.
		else
			test4(j) = .FALSE.
		end if
	end do

	call fix_shift_MN(save_i, save_m1, save_m2, save_bs, save_M, save_N, test4, cc, strict)
	if (strict) then
		 go to 9999
	end if

    do j = 2, cc
        M2(save_i(j)) = M2(save_i(j)) + save_M(j)
        N2(save_i(j)) = N2(save_i(j)) + save_N(j)
    end do 

    RETURN

9999 CONTINUE
            print*,'mm1_ori,mm2_ori, M/N, M,N | Ek, Em, En | bz, dM, dN' 
            print*, mm1_ori, mm2_ori, M/N, M, N
            print*, Ek, Em, En
            print*, (mm2_ori-mm1_ori), (sum(save_M(2:cc))-M), (sum(save_N(2:cc))-N)
            
            print*,'j, save_i, save_m1, save_m2, save_M, save_N'
            do j = 2, cc
                print*, j-1, save_i(j), save_m1(j), save_m2(j), save_bs(j), save_M(j), save_N(j), save_M(j)/save_N(j)
				print*, j-1, min_bs_ind-1, test1(j), test3(j), test4(j)
            end do
end subroutine redistribute_cloud

subroutine fix_shift_MN(save_i, save_m1, save_m2, save_bs, save_M, save_N, test4, cc,strict)
	integer, dimension(:), INTENT(IN) :: save_i
	real, dimension(:), INTENT(IN) :: save_m1, save_m2, save_bs, save_N
	real, dimension(:), INTENT(INOUT) :: save_M
	logical, dimension(:), INTENT(INOUT) :: test4
	integer, INTENT(IN) :: cc
	logical, INTENT(OUT), optional :: strict
	logical :: do_error
	real, dimension(4) :: save_M_temp
	integer :: j, j1, j2
	do_error = .FALSE.

	do j = 2, cc
		if (test4(j)) then
			if (save_M(j)/save_N(j) < save_m1(j)) then
				save_M_temp(1) = (save_m1(j) + save_bs(j) * epsL) * save_N(j)
			else
				save_M_temp(1) = (save_m2(j) - save_bs(j) * epsL) * save_N(j)
			end if

			save_M_temp(2) = save_M_temp(1) - save_M(j)
			save_M_temp(3) = save_M(j-1) - save_M_temp(2)
			save_M_temp(4) = save_M(j+1) - save_M_temp(2)

			call locate_m_new(save_M_temp(3)/save_N(j-1), j1,& 
					save_i(j-1), save_i(j-1),.FALSE.)
			call locate_m_new(save_M_temp(4)/save_N(j+1), j2,& 
					save_i(j+1), save_i(j+1),.FALSE.)

			if ((j1 == -999) .and. (j2 > -999)) then
				save_M(j+1) = save_M_temp(4)
				test4(j+1) = .FALSE.
			else if ((j1 > -999) .and. (j2 == -999)) then
				save_M(j-1) = save_M_temp(3)
				test4(j-1) = .FALSE.
			else if ((j1 == -999) .and. (j2 == -999)) then
				if (present(strict)) then
					print*, 'Warning: failed in test4 nearby shift'
					print*, 'save_M_temp',save_M_temp
					do_error = .TRUE.
					return
				else
					print*, 'Warning: failed in test4 nearby shift'
				end if
			else
				if (save_N(j+1)>save_N(j-1)) then
					save_M(j+1) = save_M_temp(4)
					test4(j+1) = .FALSE.
				else
					save_M(j-1) = save_M_temp(4)
					test4(j+1) = .FALSE.
				end if
			end if

			save_M(j) = save_M_temp(1)
			test4(j) = .FALSE.

		end if
	end do
end subroutine fix_shift_MN


subroutine line2MN(m1, m2, Ek, Em, En, dM, dN)
real, INTENT(IN) :: m1, m2, Ek, Em, En
real, INTENT(OUT) :: dM, dN
real :: dM1, dM2, dM3

dN = ( m2 - m1 ) * (En - Ek * (Em - (m1 + m2) / 2.e0))
dM1 = En * (m2**2.e0 - m1**2.e0) / 2.e0 
dM2 = (m2**3.e0 - m1**3.e0) / 3.e0
dM3 = Em * (m2**2.e0 - m1**2.e0) / 2.e0
dM = dM1 + Ek * (dM2 - dM3)

end subroutine line2MN
    
subroutine m_end_grow(mm, G, s, Af, Bf, dt,ii,jj,kk)
    real, INTENT(IN) :: G, s, dt, Af, Bf
    real, INTENT(INOUT) :: mm
	integer, INTENT(IN), optional :: ii,jj,kk
    real :: rr

    call m2r(mm,rr)

    call old_radius_liquid_euler(rr,dt,G,s,Af,Bf,ii,jj,kk)

    call r2m(rr,mm)

end subroutine m_end_grow

subroutine old_grow(rr,dt, G,s,Af,Bf,r_aero, ii,jj,kk)
	integer, INTENT(IN), optional :: ii, jj, kk
	real, intent(inout) :: rr
	real, intent(in) :: dt, G, s, Af, Bf, r_aero
	real :: dt_eff, s_eff, r_grow
	integer :: n_loop, i

    if (docs .and. (rr<5.0*rcdrop)) then
      dt_eff = 0.000001
      n_loop = nint(dt/dt_eff)
    else
      dt_eff = dt
      n_loop = 1
    end if
    do i = 1, n_loop
        if (docs) then
          s_eff = s - Af/rr + Bf/rr**3
        else
          s_eff = s
        end if
      r_grow = ((rr+r0_cor)**2.e0 + 2.e0 * G * s_eff * dt_eff)**0.5 - r0_cor
      rr = r_grow
	  if ((s<0) .and. (rr<=r_aero)) then
	  	rr = r_aero
	  	return
	  end if
    end do

	rr = max(rr, r_aero)

end subroutine old_grow

subroutine M_sum_grow(M, Mn, N, G, s, Af, Bf, dt, mm1, mm2,ii,jj,kk)
	real, INTENT(IN) :: mm1, mm2
    real, INTENT(IN) :: M, N
    real, INTENT(IN) :: G, s, dt, Af, Bf
    real, INTENT(OUT) :: Mn
	integer, INTENT(IN), optional :: ii, jj, kk
	real :: r_grow, m_grow

	if (mm1>=mm2) then
		m_grow = (mm1+mm2)/2.0e0
	else 
		r_grow = (M / N * 3.e0 / ( 4.e0 * pi * rho_H2O))**(1.e0/3.e0)
		call old_radius_liquid_euler(r_grow,dt,G,s,Af,Bf)
		call r2m(r_grow, m_grow)

		if ((mm1>m_grow) .or. (mm2<m_grow)) then
			m_grow = (mm1 + mm2)/2.0e0
		end if
	end if
    Mn = m_grow * N
end subroutine M_sum_grow

recursive subroutine old_radius_liquid_euler(r_ini,dt_int,G_pre,&
                               supersat, afactor, bfactor,ii,jj,kk)
! adapted from Fabian Hoffmann
implicit none
integer, INTENT(IN), optional :: ii, jj, kk
real, intent(inout) :: r_ini
real, intent(in) :: dt_int, G_pre, supersat, afactor, bfactor
integer :: m, cc, iii
real :: r_eul,r_eul2, r_eul_old, rel_change, dr2dt, d2r2dtdr2, f, dfdr2, r_aero 
real :: dt_eul, t_eul, r_guess

r_eul     = r_ini 
r_eul_old = r_ini

r_aero = raerosol(22) ! monodisperse aerosol

dt_eul = dt_int

cc = 0


r_guess = 1.0e-5

do m = 1, 500 

	dr2dt = 2.0 * G_pre &
		   * ( supersat - afactor / r_eul + bfactor / r_eul**3 ) &
		   * r_eul / ( r_eul + r0_cor )
	d2r2dtdr2  = G_pre & 
		   * ( afactor * r_eul**3 - bfactor * ( 3.0 * r_eul + 2.0 * r0_cor ) & 
		   + r_eul**3 * r0_cor * supersat ) & 
		   / ( r_eul**4 * ( r_eul + r0_cor )**2 )
!
!--    To speed up the Newton-Raphson scheme, the square root is executed at every iteration
	f   = r_eul**2 - r_ini**2 - dt_eul * dr2dt
	dfdr2   = 1.0 - dt_eul * d2r2dtdr2 ! = 1.0 - ( 0.0 + dt_int * d2r2dtdr2 )

	if ((r_eul**2 - f/dfdr2)<=0) then
		call old_radius_liquid_euler(r_ini,dt_int/2.0e0,G_pre,&
							   supersat, afactor, bfactor,ii,jj,kk)
		call old_radius_liquid_euler(r_ini,dt_int/2.0e0,G_pre,&
                               supersat, afactor, bfactor,ii,jj,kk)
		return
	end if

	r_eul = SQRT(r_eul**2 - f/dfdr2)

	rel_change = ABS( r_eul - r_eul_old ) / r_eul_old

	r_eul_old  = r_eul

	if (rel_change.lt.1.0e-12) then
		exit
	endif

end do
	r_ini = max(r_eul, r_aero)

	return
end subroutine old_radius_liquid_euler
		
subroutine radius_liquid_euler(r_ini,dt_int,G_pre,&
                               supersat, afactor, bfactor,ii,jj,kk)
! adapted from Fabian Hoffmann
implicit none
integer, INTENT(IN), optional :: ii, jj, kk
real, intent(inout) :: r_ini
real, intent(in) :: dt_int, G_pre, supersat, afactor, bfactor
integer :: m, cc, iii
real :: r_eul,r_eul2, r_eul_old, rel_change, dr2dt, d2r2dtdr2, f, dfdr2, r_aero 
real :: dt_eul, t_eul 

r_eul     = r_ini 
r_eul_old = r_ini

r_aero = raerosol(22) ! monodisperse aerosol

dt_eul = dt_int
t_eul = 0e0

cc = 0
iii = 0
do while (t_eul .lt. dt_int - 1.0e-20)

	do m = 1, 500 
	iii= iii+1

dr2dt = 2.0 * G_pre &
       * ( supersat - afactor / r_eul + bfactor / r_eul**3 ) &
       * r_eul / ( r_eul + r0_cor )
d2r2dtdr2  = G_pre & 
       * ( afactor * r_eul**3 - bfactor * ( 3.0 * r_eul + 2.0 * r0_cor ) & 
       + r_eul**3 * r0_cor * supersat ) & 
       / ( r_eul**4 * ( r_eul + r0_cor )**2 )

if (m==1) then
	dt_eul = min(0.5*abs(1.0/d2r2dtdr2), dt_int - t_eul, dt_int)
end if

!--    To speed up the Newton-Raphson scheme, the square root is executed at every iteration
	f   = r_eul**2 - r_ini**2 - dt_eul * dr2dt
	dfdr2   = 1.0 - dt_eul * d2r2dtdr2 ! = 1.0 - ( 0.0 + dt_int * d2r2dtdr2 )

	if ((r_eul**2 - f/dfdr2)<=0) then
		print*, 'radius: ',ii, jj, kk, dt_eul, r_eul**2 - f/dfdr2
	end if

	r_eul = SQRT( MAX( r_eul**2 - f/dfdr2, r_aero**2))

	rel_change = ABS( r_eul - r_eul_old ) / r_eul_old


	r_eul_old  = r_eul

	if (rel_change.lt.1.0e-12) then
		exit
	endif

	end do
	t_eul = t_eul + dt_eul
	r_ini = r_eul

	end do

	return
end subroutine radius_liquid_euler

subroutine get_linear_eq_E(m1, m2, M, N, Ek, Em, En)
    real, INTENT(INOUT) :: M, N
    real, INTENT(OUT) :: Ek, Em, En
    real, INTENT(INOUT) :: m1, m2
    real :: mc, bs, nc
    real :: n1n, n2n
    real :: m1o, m2o, MN_temp, M_old
    character(len=3) :: eq_name
    
    bug_happen = .FALSE.

    m1o = m1
    m2o = m2

    if (N < 0e0) then
        eq_name = 'eq0'
        if (DEBUG>0) then
            print *,'WARNING: no linear for nc1<eps'
        end if
        Ek = 0.e0
        En = 0.e0
        Em = 0.e0
    else   
        MN_temp = M/N
        if (MN_temp>m2 .or. MN_temp<m1) then 
			M_old = M
			if (MN_temp<m1) then
				M = (m1 + (m2 - m1) * epsL) * N
			else
				M = (m2 - (m2 - m1) * epsL) * N
			end if
            print*, 'Warning: get_linear_eq_E: MN',MN_temp, ' not in range m1m2'
            print*,'m1,m2',m1,m2
            print*,'M, M_old, N, newM',M,M_old, N
        end if 
                
        mc = (m1 + m2) / 2.e0
        bs = m2 - m1
        nc = N / bs
        
        Ek = 12.e0 * (M - mc * N) / bs**3.e0
        En = nc
        Em = mc

        call value_linear_eq(Ek, Em, En, m1, n1n)
        call value_linear_eq(Ek, Em, En, m2, n2n)

        if (n1n > 0e0 .and. n2n > 0e0) then
            return
        else
            En = 0.e0
            if (n1n <= 0.e0) then
                eq_name = 'eq1'
     
                Em = 3.e0 * M / N - 2.e0 * m2
                Ek = 2.e0 * N / (m2 - Em)**2.e0
                m1 = Em
            else if (n2n <= 0e0) then
                eq_name = 'eq2'
                Em = 3.e0 * M / N - 2.e0 * m1
                Ek = - 2.e0 * N / (m1 - Em)**2.e0                                          
                m2 = Em  
            else                               
                if (DEBUG>0) then
                    print *, 'Error in get_linear_eq_E'
                    print *, M, N  
                end if           
            end if
        end if
        if (m1<m1o-epsL .or. m2>m2o+epsL) then 
            print*, 'ERROR: Linear result out of original range'
            bug_happen = .TRUE.
        end if 
        if (m1>m2) then
            print*, 'ERROR: Linear result m1 m2 in the wrong order'
            bug_happen = .TRUE.
        end if 
        if (bug_happen) then
            print*, m1o, m2o, M, N
            print*, m1, m2, Ek, Em, En
            print*, eq_name, n1n, n2n
        end if
            
    end if
end subroutine get_linear_eq_E

subroutine value_linear_eq(Ek1,Em1,En1, x, y)
      real, INTENT(IN) :: Ek1, Em1, En1, x
      real, INTENT(OUT) :: y
      y = En1 + Ek1 * (x - Em1)
end subroutine value_linear_eq

subroutine calc_mass_ratio(mass_ratio)
    real :: mass_ratio
    mass_ratio = 2.e0**(1.e0/mcindex)
end subroutine calc_mass_ratio

subroutine r2m(rr,mm)
        real :: rr, mm
    mm = 4.e0 * pi / 3.e0 * rr**3 * rho_H2O
    return
end subroutine r2m

subroutine m2r(mm,rr)
        real :: rr, mm, temp
    temp = 6.e0 / pi * mm
    rr = temp**(1.e0/3.e0)/2.e0
    return
end subroutine m2r

function MN2r(M,N,kr_ind) result(rr)
    real, intent(in) :: M, N
    real :: rr
    real :: mm
    integer :: kr_ind
    if (N==0e0) then
        if (M>0e0) then
            print*,'ERROR: MN2r N=0 M>0'
        end if
        mm = (mcloud(kr_ind) + mcloud2(kr_ind))/2e0
    else
        mm = M/N
    end if
    call m2r(mm,rr)
    return
end function MN2r

function rM2N(rr, M) result(N)
    real, intent(in) :: rr, M
    real :: N
    real :: mm
    call r2m(rr,mm)
    N = M / mm
    return
end function

subroutine growth_factor(tt,qq,pp,Gcond,ss_env,A_factor,B_factor,ii,jj,kk)

  implicit none
  integer, INTENT(IN), optional :: ii, jj, kk

  real :: tt, qq, pp, Gcond, A_factor, B_factor
  real :: ss_env, esat
  real :: Dv, kt, Lv
  real :: raer_mon

  Dv = diffelq

  kt = therco

  Lv = lcond

  esat = esatw(tt)*1.e2

  ss_env = qsatw(tt,pp*1.e-3)

  ss_env = qq/ss_env

  Gcond = (rho_H2O*1.e3*rv*tt/(Dv*esat) + &
           rho_H2O*1.e3*Lv/(kt*tt) * &
           (Lv/(rv*tt)-1.0))**(-1) * 1.0e4
 
  if (docs) then
    raer_mon = raerosol(22) ! monodisperse aerosol
    A_factor = 2.0*sigma_H20/rho_H2O/rv*1.0e-4/tt
    B_factor = i_vant_Hoff(aer_type)*mol_H2O*rho_sol(aer_type) / &
                  rho_H2O / mol_weight(aer_type) * raer_mon**3
  else
    A_factor = 0e0
	B_factor = 0e0
  end if

  return
end subroutine growth_factor

subroutine aerosol_deliq2(tt,qq,pp,roro,naer,ncld,mcld,nnew,mnew)

  implicit none

  real :: tt, qq, pp, roro
  real, dimension(naerosol) :: naer, scrt, rcrt
  real, dimension(ncloud) :: ncld, mcld
  real :: nnew, mnew
  real :: qq_sat,ss_env

  call ss_environment(tt,qq,pp,qq_sat,ss_env)

  call deliq2(ss_env, tt, naer, ncld, mcld)

  return
end subroutine aerosol_deliq2

subroutine deliq2(ss_env, tt, naer, ncld, mcld)

  implicit none
  real :: ss_env, tt
  real, dimension(naerosol) :: naer
  real, dimension(ncloud) :: ncld, mcld
  real :: delRH
  integer :: ia, ic
  real :: Gf
  real :: raer_del, ra_temp
  integer :: del_index

  delRH = aer_delRH(aer_type)
  if ( ss_env > delRH ) then
  
    do ia=1, naerosol
      if (naer(ia) > 0e0) then
	    del_index = 1
		ra_temp = raerosol(ia)
		raer_del = r_equi(ss_env, tt, ra_temp)
		do ic = 1, ncloud
		  if (ic==1.and.raer_del<rcloud(ic)) then
            ncld(ic) = ncld(ic) + naer(ia)
	        mcld(ic) = mcld(ic) + 4.0*pi/3.0*rho_H2O* &
			      ((rcloud(ic)+rcloud(ic+1))/2.0)**3*naer(ia)
            naer(ia) = 0.0
			del_index = 0
	      elseif (ic>1.and.raer_del<rcloud(ic)) then
            ncld(ic-1) = ncld(ic-1) + naer(ia)
	        mcld(ic-1) = mcld(ic-1) + 4.0*pi/3.0*rho_H2O* &
			      ((rcloud(ic-1)+rcloud(ic))/2.0)**3*naer(ia)
            naer(ia) = 0.0
			del_index = 0
		  end if
		  if (del_index==0) exit
		end do
	  end if
	end do

  end if

  return
end subroutine deliq2

real function r_equi(ss,tt,r_aero)
  implicit none
  real, intent(in) :: ss, tt, r_aero
  real :: ss_temp
  real :: A_factor, B_factor, r_equi_0

    A_factor = 2.0*sigma_H20/rho_H2O/rv*1.0e-4/tt
    B_factor = i_vant_Hoff(aer_type)*mol_H2O*rho_sol(aer_type) / &
                  rho_H2O / mol_weight(aer_type)
	r_equi_0 = 1.0
	r_equi   = 1.0E-4

	ss_temp  = max(-0.2, min(ss-1.0,-0.01))

	do while (abs((r_equi - r_equi_0)/r_equi_0).gt.1.0e-20)
	  r_equi_0 = r_equi
	  r_equi   = ((B_factor*r_aero**3)/(A_factor/r_equi_0-ss_temp))**(1.0/3.0)
	end do
  
end function r_equi

subroutine aerosol_deliq(tt,qq,pp,roro,naer,ncld,mcld,nnew,mnew)

  implicit none

  real :: tt, qq, pp, roro
  real, dimension(naerosol) :: naer, scrt, rcrt
  real, dimension(ncloud) :: ncld, mcld
  real :: nnew, mnew
  real :: qq_sat,ss_env

  call ss_environment(tt,qq,pp,qq_sat,ss_env)

  call deliq(ss_env, naer, ncld, mcld)

  return
end subroutine aerosol_deliq

subroutine deliq(ss_env, naer, ncld, mcld)

  implicit none
  real :: ss_env
  real, dimension(naerosol) :: naer
  real, dimension(ncloud) :: ncld, mcld
  real :: delRH
  integer :: ia, ic
  real :: Gf
  real :: raer_del
  integer :: del_index

  delRH = aer_delRH(aer_type)
  if ( ss_env > delRH ) then
  
    do ia=1, naerosol
      if (naer(ia) > 0e0) then
	    del_index = 1
	    call G_factor(ia, Gf)
		raer_del = Gf * raerosol(ia)
		do ic = 1, ncloud
		  if (ic==1.and.raer_del<rcloud(ic)) then
            ncld(ic) = ncld(ic) + naer(ia)
	        mcld(ic) = mcld(ic) + 4.0*pi/3.0*rho_H2O* &
			      ((rcloud(ic)+rcloud(ic+1))/2.0)**3*naer(ia)
            naer(ia) = 0.0
			del_index = 0
	      elseif (ic>1.and.raer_del<rcloud(ic)) then
            ncld(ic-1) = ncld(ic-1) + naer(ia)
	        mcld(ic-1) = mcld(ic-1) + 4.0*pi/3.0*rho_H2O* &
			      ((rcloud(ic-1)+rcloud(ic))/2.0)**3*naer(ia)
            naer(ia) = 0.0
			del_index = 0
		  end if
		  if (del_index==0) exit
		end do
	  end if
	end do

  end if

  return
end subroutine deliq

subroutine G_factor(ia, Gf)

  implicit none
  integer :: ia
  real :: Gf
  real :: delRH, ivff, rhoa, Maer, Mwater

  delRH = aer_delRH(aer_type)
  ivff = i_vant_Hoff(aer_type)
  rhoa = rho_sol(aer_type)
  Maer = mol_weight(aer_type)
  Mwater = mol_H2O

  Gf = ((1.0-delRH+(ivff*rhoa*Mwater/rho_H2O/Maer)*delRH) &
           / (1.0-delRH))**(1d0/3d0)
  return
end subroutine G_factor

subroutine activation(ss_env, scrt, rcrt, naer, ncld, mcld, nnew, mnew)

  implicit none

  real :: ss_env
  real,dimension(naerosol) :: naer, scrt, rcrt
  real,dimension(ncloud) :: ncld, mcld
  real :: nnew, mnew
  integer :: ia

  do ia = naerosol,1,-1
    if (scrt(ia)<ss_env) then
      ncld(1) = ncld(1) + naer(ia)
      mcld(1) = mcld(1) + 4.0*pi/3.0*rho_H2O* &
                ((rcloud(1)+rcloud(2))/2.0)**3*naer(ia)
      nnew = nnew + naer(ia)
      mnew = mnew + 4.0*pi/3.0*rho_H2O* &
             ((rcloud(1)+rcloud(2))/2.0)**3*naer(ia)
      naer(ia) = 0.0
    else
      exit
    end if
  end do
  return
end subroutine activation


subroutine activation2(tt, roro, qq, pp, ss_env, scrt, rcrt, naer, ncld, mcld, nnew, mnew)

  implicit none

  real :: tt, roro, qq, pp, ss_env
  real,dimension(naerosol) :: naer, scrt, rcrt
  real,dimension(ncloud) :: ncld, mcld
  real :: nnew, mnew
  integer :: ia, i

  real :: maxql, naer2cld, maxnact, ttemp, qqtemp, qq_sat, sstemp

  ttemp = tt
  qqtemp = qq
  nnew = 0.0
  mnew = 0.0
  qq_sat = qsatw(ttemp,pp*1.e-3)
  sstemp = qqtemp/qq_sat
  do ia = naerosol,1,-1
    if (scrt(ia)<sstemp) then
      maxql = qsatw(ttemp,pp*1.e-3)*(ss_env - scrt(ia))
      maxnact = maxql*roro/m1dropnew - nnew
      if (maxnact < 0.0) then
        exit
      elseif (0.2*maxnact>naer(ia)) then
        naer2cld = naer(ia)
      else
        naer2cld = maxnact*0.2 ! correction factor
      end if
      ncld(1) = ncld(1) + naer2cld
      mcld(1) = mcld(1) + m1dropnew*naer2cld
      nnew = nnew + naer2cld
      mnew = mnew + m1dropnew*naer2cld
      naer(ia) = naer(ia) - naer2cld
      ttemp = ttemp + fac_cond*m1dropnew*naer2cld/roro
      qqtemp = qqtemp - m1dropnew*naer2cld/roro
      qq_sat = qsatw(ttemp,pp*1.e-3)
      sstemp = qqtemp/qq_sat
      if (naer(ia)>1.0e-15) then
        exit
      end if
    else
      exit
    end if
  end do
  return
end subroutine activation2

subroutine ss_environment(tt,qq,pp,qq_sat,ss_env)

  implicit none

  real :: tt, qq, pp
  real :: qq_sat, ss_env

  qq_sat = qsatw(tt,pp*1.e-3) 
 
  ss_env = qq/qq_sat

  return

end subroutine ss_environment

subroutine ss_critical(tt,scrt,rcrt)

  implicit none

  real :: tt
  real, dimension(naerosol) :: raer, scrt, rcrt

  integer :: ia
  real :: scritical_A, scritical_B
  real :: rincm

  if (aer_type == 1) then
    scritical_A = 2.0*sigma_H20/rho_H2O/rv*1.0e-4/tt
    scritical_B = i_vant_Hoff(aer_type)*mol_H2O*rho_sol(aer_type) / &
                  rho_H2O / mol_weight(aer_type)
  elseif (aer_type == 2) then
    scritical_A = 2.0*sigma_H20/rho_H2O/rv*1.0e-4/tt
    scritical_B = i_vant_Hoff(aer_type)*mol_H2O*rho_sol(aer_type) / &
                  rho_H2O / mol_weight(aer_type)
  else
    print*,"Error: Aerosol Type Not Defined"
    stop
  end if

  do ia = 1, naerosol
    rincm = raerosol(ia)
    scrt(ia) = 1.0 + sqrt(4.0*scritical_A**3/27.0/scritical_B/rincm**3)
    rcrt(ia) = sqrt(3.0*scritical_B*rincm**3/scritical_A)
  end do

  return
end subroutine ss_critical



subroutine test
implicit none
print*,"test"
end subroutine test

subroutine get_new_bin_end(Nnew,bin_end)
real, dimension(ncloud), INTENT(IN) :: Nnew
integer, INTENT(INOUT) :: bin_end
integer :: i

    if (bin_end .eq. ncloud) then
        RETURN
    else
        do i = bin_end+1,ncloud+1
            if (i .gt. ncloud) then
                bin_end = ncloud
                RETURN
            else
                if (Nnew(i) .lt. eps60) then
                    bin_end = i-1
                    RETURN
                end if
            end if
        end do
    end if
end subroutine get_new_bin_end
        
subroutine get_new_bin_begin(Nnew,bin_begin)
real, dimension(ncloud), INTENT(IN) :: Nnew
integer, INTENT(INOUT) :: bin_begin
integer :: i
    if (bin_begin .eq. 1) then
        RETURN
    else
        do i = bin_begin-1,0,-1
            if (i .lt. 1) then
                bin_begin = 1
                RETURN
            else
                if (Nnew(i) .lt. eps60) then
                    bin_begin = i+1
                    RETURN
                end if
            end if
        end do
    end if
end subroutine get_new_bin_begin

end module module_ntubm
