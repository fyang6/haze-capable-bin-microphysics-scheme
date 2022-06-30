module module_ntubm

! This module is a double-moment bin microphysical schemes based on Chen and
! Lamb JAS 1994
! Reference:
! Chen, J.P. and Lamb, D., 1994. Simulation of cloud microphysical and chemical
! processes using a multicomponent framework. Part I: Description of the
! microphysical model. Journal of the Atmospheric Sciences, 51(18),
! pp.2613-2630.
!
! The code is implemented in SAM by Fan Yang @ BNL, 2021
!
! History:
! version, date, description
! v0.1(F. Yang); 01/12/2021; only warm cloud microphysical processes are
! included
! v0.2(F. Yang); 06/21/2021; bug fixed

use params, only: cp, ggr, rgas, rv, fac_cond, lcond, pi, diffelq, therco 
use grid, only: masterproc, rank
use vars, only: esatw, qsatw
use micro_prm

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
        
subroutine cloud_collision(dtcond, ncld, mcld, bin_begin, bin_end)
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

    if (abs(sum(deltaM))>eps) then
        print*, 'ERROR: after collision, M not balanced'
        print*, 'deltaM:', deltaM
        call sleep(2)
        stop
    end if

    if (sum(deltaN)>0) then
        print*, 'ERROR: N increase after collision'
        call sleep(2)
        stop
    end if
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

subroutine cloud_cond_growth(tt, qq, pp, roro, dtcond, ncld, mcld, ncnew, bin_begin, bin_end) 
    real, INTENT(IN) :: pp, roro, dtcond
    integer, INTENT(INOUT) :: bin_begin, bin_end
    real, INTENT(INOUT) :: tt, qq
    real, INTENT(OUT) :: ncnew
    real, dimension(ncloud), INTENT(INOUT) :: ncld, mcld
    real, dimension(ncloud) :: Nnew, Mnew
    real :: Gcond, ss_env
    real, dimension(ncloud) :: N_temp
    real, dimension(ncloud) :: M_N_temp
    real :: M_Ni, m1_temp,m2_temp, Em, Ek, En
    integer :: i
    
    ncnew = 0.0

    Nnew = 0.0
    Mnew = 0.0

    call growth_factor(tt,qq,pp,Gcond,ss_env)  
    do i = bin_begin, bin_end
        if (ncld(i) > 0 ) then
            call do_one_cloud_growth(i, mcld(i), ncld(i), Mnew, Nnew, ncnew, Gcond, ss_env, dtcond)
        end if 
    end do

    call get_new_bin_begin(Nnew,bin_begin)
    call get_new_bin_end(Nnew,bin_end)

    ncld = Nnew
    mcld = Mnew
end subroutine cloud_cond_growth

subroutine do_one_cloud_growth(i, M1, N1, M2, N2, N3, G, s_1, dt)
    integer, INTENT(IN) :: i
    real, INTENT(IN) :: G, s_1, dt
    real, INTENT(INOUT) :: M1, N1
    real, dimension(ncloud), INTENT(INOUT):: M2, N2
    real, INTENT(INOUT) :: N3
    real :: mm1, mm2, Mn, Ek, Em, En
    real :: ss
    logical :: debug = .False.

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
    
    call m_end_grow(mm1, G, ss, dt) 
    call m_end_grow(mm2, G, ss, dt) 
    if (mm1>=mcloud(1)) then
       call M_sum_grow(M1, Mn, N1, G, ss, dt) 

       call get_linear_eq_E(mm1, mm2, Mn, N1, Ek, Em, En) 
       if (bug_happen) then
          print*, 'bug_happen: after 2nd linear in do_one'
          call sleep(2)
          stop
       end if
       call redistribute_cloud(mm1, mm2, Mn, N1, Ek, Em, En, M2, N2, i)
    else
       N3 = N3 + N1  
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
    integer :: left_i, right_i

    left_i = guess
    right_i = guess

    do while (mm1_ori < mcloud(left_i))
        left_i = left_i - 1
    end do
    do while (mm2_ori > mcloud2(right_i))
        right_i = right_i + 1
    end do
    do while (mm1_ori > mcloud2(left_i))
        left_i = left_i + 1
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
        print*, 'Error: redistribution find wrong ind edge', left_i, right_i, mm1_ori, mm2_ori
        call sleep(2)
        stop
    end if

    ! Case 2: exceed edge, then ignore all 
    if (left_i<1) then
        print*, 'WARNING: evaportae', mm1_ori, mm2_ori, M, N
        RETURN
    else if (right_i>ncloud) then
        print *, 'WARNING: too large',mm1_ori, mm2_ori, M, N
        RETURN
    end if

    ! Case 3: extreme small bin, or very few N
    if (mm2_ori-mm1_ori<eps .or. N<eps) then
        call locate_m_new(M/N, j2, left_i, right_i,.true.)
        M2(j2) = M2(j2) + M
        N2(j2) = N2(j2) + N
        RETURN
    end if 

    ! Case 4: distribute to left_i ~ right_i bins
    cc = 0
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
        call line2MN(save_m1(cc), save_m2(cc), Ek, Em, En, save_M(cc), save_N(cc))
        save_bs(cc) = save_m2(cc) - save_m1(cc)
        save_i(cc) = i
    end do
        
    ! Multiple tests for Case 4
    dM_temp = sum(save_M(1:cc))
    dN_temp = sum(save_N(1:cc))
    if ((dM_temp < M * (1e0 - epsR)) .or. &
        (dM_temp > M * (1e0 + epsR)) .or. &
        (dN_temp < N * (1e0 - epsR)) .or. &
        (dN_temp > N * (1e0 + epsR))) then

        if ((mm2_ori-mm1_ori<eps*10e0) .or. &
            (abs(dN_temp-N)<eps*1e-5) .or. &
            (abs(Ek)>1e10) .or. &
            (abs(dN_temp-N)<eps .and. abs(Ek)>1e5)) then

            if (cc>2) then
                print*, 'Error: redistribution Test 2: extreme small bin cannot seperate more than 2 fixed bin'
                go to 9999
            end if 

            call locate_m_new(M/N, j2, left_i, right_i, .true.)
            M2(j2) = M2(j2) + M
            N2(j2) = N2(j2) + N
            RETURN
        else
            if (save_bs(1) < save_bs(cc)) then
                if (save_bs(1)<eps) then
                    save_M(1) = M - sum(save_M(2:cc))
                    save_N(1) = N - sum(save_N(2:cc))
                else
                    print*, 'ERROR: redistribute Test1.2a more or less M or N'
                    go to 9999
                end if
            else
                if (save_bs(cc)<eps) then
                    save_M(cc) = M - sum(save_M(1:(cc-1)))
                    save_N(cc) = N - sum(save_N(1:(cc-1)))
                else
                    print*, 'ERROR: redistribute Test1.2b more or less M or N'
                    go to 9999
                end if
            end if 
        end if
    end if

    ! Test 3: test negative value for both edge
    if (save_M(1)<0 .or. save_N(1)<0) then
        save_M(1) = 0
        save_N(1) = 0
        save_M(cc) = M - sum(save_M(2:(cc-1)))
        save_N(cc) = N - sum(save_N(2:(cc-1)))
    end if

    if (save_M(cc)<0 .or. save_N(cc)<0) then
        save_M(cc) = 0
        save_N(cc) = 0
        save_M(1) = M - sum(save_M(2:cc))
        save_N(1) = N - sum(save_N(2:cc))
    end if
            
    ! Test 4: test M/N locate
    do j = 1, cc
       call locate_m_new(save_M(j)/save_N(j), j1, save_i(j), save_i(j),.FALSE.)
       if (j1 == -999) then 
            if (save_N(j)<eps) then
                if (j==1) then
                    save_M(1) = 0
                    save_N(1) = 0
                    if (c==2) then
                        save_M(2) = M
                        save_N(2) = N
                    else if (c>2) then
                        save_M(2) = M - sum(save_M(3:cc))
                        save_N(2) = N - sum(save_N(3:cc))
                    else! c==1
                        call locate_m_new(M/N, j1, left_i, right_i, .TRUE.)
                        M2(j1) = M2(j1) + M
                        N2(j1) = N2(j1) + N
                        RETURN
                    end if
                else if (j==cc) then
                    save_M(cc) = 0
                    save_N(cc) = 0
                    if (c==2) then
                        call locate_m_new(M/N, j1, save_i(1), save_i(1), .FALSE.)
                        if (j1==-999) then
                            if (N<eps) then
                                call locate_m_new(M/N, j1, left_i, right_i, .TRUE.)
                                if (j1==-999) then
                                    print*, 'Error: redistribution Test 4.1a cannot relocate shifted eps value'
                                    go to 9999
                                else
                                    M2(j1) = M2(j1) + M
                                    N2(j1) = N2(j1) + N
                                    RETURN
                                end if
                            else
                                print*, 'Error: redistribution Test 4.1b nearby bin shifted with eps added'
                                go to 9999
                            end if
                        end if
                    else if (c>2) then
                        save_M(cc-1) = M - sum(save_M(1:(cc-2)))
                        save_N(cc-1) = N - sum(save_N(1:(cc-2)))
                        call locate_m_new(save_M(cc-1)/save_N(cc-1), j1, save_i(cc-1), save_i(cc-1), .FALSE.)
                        if (j1==-999) then
                            if (save_N(cc-1)<eps) then
                                call locate_m_new(save_M(cc-1)/save_N(cc-1), j1, &
                                        left_i, right_i, .TRUE.)
                                if (j1==-999) then
                                    print*, 'Error: redistribution Test 4.1c cannot relocate shifted eps value'
                                    go to 9999
                                else
                                    save_i(cc-1) = j1
                                end if
                            else
                                print*, 'Error: redistribution Test 4.1d nearby bin shifted with eps added'
                                go to 9999
                            end if
                        end if
                    end if
                else
                    print*, 'Error: redistribution Test 4.1e: eps in the middle'
                    go to 9999
                end if
            else ! Test 4.2: if not neglegible
                m1n = save_M(j)/save_N(j)
                m2n = mcloud2(save_i(j))-mcloud(save_i(j)) 
                if (mcloud(save_i(j))>m1n) then 
                    div_temp = (mcloud(save_i(j))-m1n)/m2n*100e0 
                else ! shift to right
                    div_temp = (m1n-mcloud2(save_i(j)))/m2n*100e0
                end if
                if (div_temp > 10.0e0) then
                    print*, 'WARNING: redistribute Test 4.3 shift',div_temp, '% of bin size', m2n
                    print*, 'm1,m2',mcloud(save_i(j)), mcloud2(save_i(j))
                    print*, 'M, N, M/N',save_M(j), save_N(j), m1n
                    print*, 'Em, En, Ek', Em, En, Ek
                end if
                if (div_temp>50e0) then
                    if (abs(Ek)>1e8) then
                        call locate_m_new(M/N, j2, left_i, right_i, .TRUE.)
                        M2(j2) = M2(j2) + M
                        N2(j2) = N2(j2) + N
                        RETURN
                    else
                        print*, 'ERROR: Test 4.3 shift too much'
                        go to 9999
                    end if
                end if
                call locate_m_new(save_M(j)/save_N(j), j1, left_i, right_i,.FALSE.)
                save_i(j) = j1
            end if
        end if
    end do

    do j = 1, cc
        M2(save_i(j)) = M2(save_i(j)) + save_M(j)
        N2(save_i(j)) = N2(save_i(j)) + save_N(j)
    end do 

    RETURN

9999    CONTINUE
            print*,'mm1_ori,mm2_ori,M,N | Ek, Em, En | bz, dM, dN' 
            print*, mm1_ori, mm2_ori, M, N
            print*, Ek, Em, En
            print*, (mm2_ori-mm1_ori), (dM_temp-M), (dN_temp-N)
            
            print*,'j, save_i, save_m1, save_m2, save_M, save_N'
            do j = 1, cc
                print*, j, save_i(j), save_m1(j), save_m2(j), save_M(j), save_N(j)
            end do
            call sleep(2)
            stop
end subroutine redistribute_cloud


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
    
subroutine m_end_grow(mm, G, s, dt)
    real, INTENT(IN) :: G, s, dt
    real, INTENT(INOUT) :: mm
    real :: rr, r_grow
    
    call m2r(mm,rr)
    r_grow = (rr**2.e0 + 2.e0 * G * s * dt)**0.5
    call r2m(r_grow,mm)

end subroutine m_end_grow
    
subroutine M_sum_grow(M, Mn, N, G, s, dt)
    real, INTENT(IN) :: M
    real, INTENT(IN) :: N, G, s, dt
    real, INTENT(OUT) :: Mn
    real :: r_star, r_grow, m_grow
    r_star = (M / N * 3.e0 / ( 4.e0 * pi * rho_H2O))**(1.e0/3.e0)
    r_grow = (r_star**2.e0 + 2.e0 * G * s * dt)**0.5
    call r2m(r_grow, m_grow)
    Mn = m_grow * N
end subroutine M_sum_grow

subroutine get_linear_eq_E(m1, m2, M, N, Ek, Em, En)
    real, INTENT(IN) :: M, N
    real, INTENT(OUT) :: Ek, Em, En
    real, INTENT(INOUT) :: m1, m2
    real :: mc, bs, nc
    real :: n1n, n2n
    real :: m1o, m2o, MN_temp
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
            print*, 'get_linear_eq_E--ERROR: MN',MN_temp, ' not in range m1m2'
            print*,'m1,m2',m1,m2
            print*,'M, N',M,N
            bug_happen = .TRUE.
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
        if (m1<m1o-eps .or. m2>m2o+eps) then 
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

subroutine r2m(rr,mm)!L
        real :: rr, mm
    mm = 4.e0 * pi / 3.e0 * rr**3 * rho_H2O
    return
end subroutine r2m


subroutine m2r(mm,rr)!L
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

subroutine growth_factor(tt,qq,pp,Gcond,ss_env)

  implicit none

  real :: tt, qq, pp, Gcond
  real :: ss_env, esat
  real :: Dv, kt, Lv

  Dv = diffelq

  kt = therco

  Lv = lcond

  esat = esatw(tt)*1.e2

  ss_env = qsatw(tt,pp*1.e-3)

  ss_env = qq/ss_env

  Gcond = (rho_H2O*1.e3*rv*tt/(Dv*esat) + &
           rho_H2O*1.e3*Lv/(kt*tt) * &
           (Lv/(rv*tt)-1.0)*ss_env)**(-1) * 1.0e4
  
  return
end subroutine growth_factor

subroutine cloud_activation(tt,qq,pp,roro,naer,ncld,mcld,nnew,mnew)

  implicit none

  real :: tt, qq, pp, roro
  real, dimension(naerosol) :: naer, scrt, rcrt
  real, dimension(ncloud) :: ncld, mcld
  real :: nnew, mnew
  real :: qq_sat,ss_env

  call ss_environment(tt,qq,pp,qq_sat,ss_env)

  call ss_critical(tt,scrt,rcrt)

  call activation2(tt, roro, qq, pp, ss_env, scrt, rcrt, naer, ncld, mcld, nnew, mnew)
  !call activation(ss_env, scrt, rcrt, naer, ncld, mcld, nnew, mnew)
 
  return
end subroutine cloud_activation

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
  real :: A_factor, B_factor
  real :: qq_sat, ss_env, e_sat

! for HUJI comparison
  A_factor = 2.53e12
  B_factor = 5.42e3

  e_sat = A_factor * exp(-B_factor/tt)

  ss_env = qq*pp/(0.622+0.378*qq)

  ss_env = ss_env/e_sat
  
  !if (masterproc) print*,"test",esatw(tt),qsatw(tt,pp*1.e-3),e_sat,ss_env

  !qq_sat = qsatw(tt,pp*1.e-3) ! consistent 
 
  !ss_env = qq/qq_sat

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
