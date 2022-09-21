module bulk_sub
  use param
  ! -------------------------------------------------------------
  ! Modules for calculating heat,freshwater, and momentum fluxes
  ! -------------------------------------------------------------
  implicit none
  real(idx),parameter :: lh_to_ev=-4.1e-10_idx
  private
  public :: diurnal_sw
  public :: cal_ws
  public :: ld_to_lw
  public :: bulk_kara00,bulk_kara05,bulk_ncep
  public :: lh_to_ev
contains
  !==================================================
  ! Calculate wind stress from wind speeds
  ! (with constant bulk coefficient)
  !==================================================
  subroutine cal_ws(u10,v10,ws10,tau_x,tau_y)
    implicit none
    real(idx),intent(in) :: u10,v10,ws10
    real(idx),intent(inout) :: tau_x,tau_y
    real(idx) :: rho_a,C_d
    ! Setting parameters
    C_d = 1.43e-3_idx
    rho_a = 1.225_idx ! [kg/m**3]
    tau_x= rho_a * C_d * ws10 * u10
    tau_y= rho_a * C_d * ws10 * v10
  end subroutine cal_ws
  !==================================================
  ! Calculate turbulent fluxes from atmospheric data using
  ! bulk formulae of Kara et al. (2000)
  ! " Efficient and Accurate Bulk Parameterizations of
  ! Air–Sea Fluxes for Use in General Circulation Models"
  !==================================================
  subroutine bulk_kara00(t_a,q_a,t_s,ws10_in,sh,lh,ev)
    implicit none
    real(idx),intent(in) :: t_a,q_a,t_s,ws10_in
    real(idx),intent(inout) :: sh,lh,ev
    real(idx) :: ws10
    real(idx) :: e_s,q_s,L
    real(idx) :: rho_a,p_air,R_g
    real(idx) :: C_e,C_h,C_0,C_1,C_pa
    real(idx) :: t_diff
    ! Setting parameters
    p_air = 1.013e3_idx ! hPa
    R_g = 2.871e2_idx  ! J kg^(-1) K^(-1)
    C_pa = 1.00467e3_idx  !  J kg^(-1) K (^-1)
    e_s=0.98_idx*((1.0007_idx+3.46e-6_idx*p_air)*6.1121_idx*exp(17.502_idx*t_s/(240.97+t_s)))
    q_s = 0.62197_idx * e_s / (p_air - 0.378_idx * e_s) ! [g/g]
    L = (1.0e6_idx) * (2.501 - 0.00237_idx * t_s) ! [J/kg]
    rho_a = (1.0e2_idx) * p_air / (R_g * (T_a + 273.16_idx) * (1.0_idx + 0.61_idx * q_a)) ! [kg/m^3]
    if (ws10_in <= 2.5_idx) then
       ws10 = 2.5_idx
    else if (ws10_in >= 32.5_idx) then
       ws10 = 32.5_idx
    else
       ws10 = ws10_in
    end if
    t_diff = t_a-t_s-(t_a+273.16_idx)*0.61_idx*(q_s-q_a)    
    C_0 = (1.0e-3_idx) * (9.94e-1_idx + (6.1e-2_idx) * ws10 - (1.0e-3_idx)*(ws10**2))
    C_1 = (1.0e-3_idx) * (-2.0e-2_idx + 6.91e-1_idx / ws10 - (8.17e-1_idx) / (ws10**2))
    C_e = C_0 + C_1 * t_diff
    C_H = 0.96_idx * C_e
    sh = -1.0_idx * rho_a * C_pa * C_H * ws10 * (T_s -T_a)
    lh = -1.0_idx * L * rho_a * C_e * ws10 * (q_s-q_a)
    ev = lh_to_ev * lh
!    ev = rho_a * C_e * ws10 * (q_s-q_a)*1.0e-3_idx
  end subroutine bulk_kara00
  subroutine bulk_kara05(t_a,q_a,t_s,ws10,sh,lh,ev)
    implicit none
    real(idx),intent(in) :: t_a,q_a,t_s,ws10
    real(idx),intent(inout) :: sh,lh,ev
    real(idx) :: t_diff
    real(idx) :: C_L0,C_L1,C_L2
    real(idx) :: C_D0,C_D1,C_D2
    real(idx) :: cmat_cd_1(3,4),cmat_cd_2(2,4),cmat_cd_3(3,4)
    real(idx) :: cmat_cd_4(3,4),cmat_cd_5(2,4),cmat_cd_6(3,4)
    real(idx) :: cmat_ch_1(3,3),cmat_ch_2(2,3),cmat_ch_3(3,3)
    real(idx) :: cmat_ch_4(3,3),cmat_ch_5(2,5),cmat_ch_6(3,3)
    real(idx) :: rho_a,p_air,R_g,C_pa,e_s,q_s,L
    real(idx) :: C_e,C_h
    p_air = 1.013e3_idx ! hPa
    R_g = 2.871e2_idx  ! J kg^(-1) K^(-1)
    C_pa = 1.00467e3_idx  !  J kg^(-1) K (^-1)
    e_s=0.98_idx*((1.0007_idx+3.46e-6_idx*p_air)*6.1121_idx*exp(17.502_idx*t_s/(240.97+t_s)))
    q_s = 0.62197_idx * e_s / (p_air - 0.378_idx * e_s) ! [g/g]
    L = (1.0e6_idx) * (2.501 - 0.00237_idx * t_s) ! [J/kg]
    rho_a = (1.0e2_idx) * p_air / (R_g * (T_a + 273.16_idx) * (1.0_idx + 0.61_idx * q_a)) ! [kg/m^3]

    ! Mat 1
    cmat_cd_1(1,1) =  1.891_idx ; cmat_cd_1(1,2) = -7.182e-1_idx
    cmat_cd_1(1,3) =  1.975e-1 ; cmat_cd_1(1,4) = -1.79e-2_idx
    cmat_cd_1(2,1) = -6.3e-3_idx ; cmat_cd_1(2,2) = -3.028e-1_idx
    cmat_cd_1(2,3) = 3.12e-1_idx ; cmat_cd_1(2,4) = -1.21e-1_idx
    cmat_cd_1(3,1) = 4.4e-4_idx ; cmat_cd_1(3,2) = -1.769e-1_idx
    cmat_cd_1(3,3) = 1.303e-2_idx ; cmat_cd_1(3,4) = -3.39e-3_idx    
    cmat_ch_1(1,1) =  2.077_idx ; cmat_ch_1(1,2) = -3.933e-1_idx ; cmat_ch_1(1,3) =  3.971e-2_idx
    cmat_ch_1(2,1) = -2.899e-1_idx ; cmat_ch_1(2,2) = 7.35e-2_idx  ; cmat_ch_1(2,3) = -6.27e-3_idx
    cmat_ch_1(3,1) = -1.954e-2_idx ; cmat_ch_1(3,2) = 5.483e-3_idx  ; cmat_ch_1(3,3) = -4.867e-4_idx

    ! Mat 2
    cmat_cd_2(1,1) =  9.774e-1_idx ; cmat_cd_2(1,2) = -2.566e-1_idx
    cmat_cd_2(1,3) =  1.048e-1_idx ; cmat_cd_2(1,4) =-1.097e-2_idx
    cmat_cd_2(2,1) = 2.051e-1_idx ; cmat_cd_2(2,2) = -1.903_idx
    cmat_cd_2(2,3) = 1.133_idx ; cmat_cd_2(2,4) = -2.658e-1_idx
    cmat_ch_2(1,1) =  8.58e-1_idx ; cmat_ch_2(1,2) = 9.743e-2_idx ; cmat_ch_2(1,3) = -1.056e-2_idx
    cmat_ch_2(2,1) = -1.927_idx ; cmat_ch_2(2,2) = 7.345e-1_idx   ; cmat_ch_2(2,3) = -7.706e-2_idx
    ! Mat 3
    cmat_cd_3(1,1) = -6.695e-2_idx ; cmat_cd_3(1,2) = 3.133e-1_idx
    cmat_cd_3(1,3) =  1.47e-3_idx ; cmat_cd_3(1,4) = -4.06e-3_idx
    cmat_cd_3(2,1) = 9.966-2_idx ; cmat_cd_3(2,2) = -2.116_idx
    cmat_cd_3(2,3) = 4.626_idx ; cmat_cd_3(2,4) = -2.68_idx
    cmat_cd_3(3,1) = -2.477e-2_idx ; cmat_cd_3(3,2) = 2.726e-1_idx
    cmat_cd_3(3,3) = -5.558e-1_idx ; cmat_cd_3(3,4) = 3.1390e-1_idx
    cmat_ch_3(1,1) = -2.925e-1_idx ; cmat_ch_3(1,2) =  5.498e-1_idx ; cmat_ch_3(1,3) = -5.544e-2_idx
    cmat_ch_3(2,1) =  7.372e-2_idx ; cmat_ch_3(2,2) = -1.74e-1_idx  ; cmat_ch_3(2,3) = 2.489e-2_idx
    cmat_ch_3(3,1) = -6.95e-3_idx ; cmat_ch_3(3,2) =  1.637e-2_idx ; cmat_ch_3(3,3) = -2.62e-3_idx
    ! Mat 4
    cmat_cd_4(1,1) =  6.497e-1_idx ; cmat_cd_4(1,2) = 6.693e-2_idx
    cmat_cd_4(1,3) =  3.54e-5_idx ; cmat_cd_4(1,4) = -3.43e-6_idx
    cmat_cd_4(2,1) = 3.83e-3_idx ; cmat_cd_4(2,2) = -2.756e-1_idx
    cmat_cd_4(2,3) = -1.091_idx ; cmat_cd_4(2,4) = 4.946_idx
    cmat_cd_4(3,1) = -4.83e-5_idx ; cmat_cd_4(3,2) = 7.71e-3_idx
    cmat_cd_4(3,3) = -2.555e-1_idx ; cmat_cd_4(3,4) = 7.654e-1_idx
    cmat_ch_4(1,1) =  1.074_idx    ; cmat_ch_4(1,2) =  5.58e-3_idx ; cmat_ch_4(1,3) =  5.26e-5_idx
    cmat_ch_4(2,1) =  6.91e-3_idx ; cmat_ch_4(2,2) = -2.244e-1_idx ; cmat_ch_4(2,3) = -1.027_idx
    cmat_ch_4(3,1) =  1.9e-4_idx ; cmat_ch_4(3,2) = -2.18e-3_idx ; cmat_ch_4(3,3) = -1.01e-1_idx
    ! Mat 5
    cmat_cd_5(1,1) =  5.438e-1_idx ; cmat_cd_5(1,2) = 8.316e-2_idx
    cmat_cd_5(1,3) =  -4.9e-4_idx ; cmat_cd_5(1,4) = 3.09e-6_idx
    cmat_cd_5(2,1) = 1.669e-2_idx ; cmat_cd_5(2,2) = 5.738e-1_idx
    cmat_cd_5(2,3) = -1.224e1_idx ; cmat_cd_5(2,4) = 3.253e1_idx
    cmat_ch_5(1,1) =  1.023395_idx ; cmat_ch_5(1,2) = 9.61e-3_idx  ; cmat_ch_5(1,3) =  -2.16e-5_idx
    cmat_ch_5(2,1) = -3.93e-3_idx  ; cmat_ch_5(2,2) = 2.048e-1_idx ; cmat_ch_5(2,3) = -5.048_idx
    ! Mat 6
    cmat_cd_6(1,1) = 5.581e-1_idx ; cmat_cd_6(1,2) = 8.174e-2_idx
    cmat_cd_6(1,3) =  -4.5e-4_idx ; cmat_cd_6(1,4) = 2.67e-6_idx
    cmat_cd_6(2,1) = -5.59e-3_idx ; cmat_cd_6(2,2) = 2.096e-1_idx
    cmat_cd_6(2,3) = -8.634_idx ; cmat_cd_6(2,4) = 1.863e1_idx
    cmat_cd_6(3,1) = 6.0e-4_idx ; cmat_cd_6(3,2) = -2.629e-2_idx
    cmat_cd_6(3,3) = 2.121e-1_idx ; cmat_cd_6(3,4) = 7.755e-1_idx
    cmat_ch_6(1,1) = 1.023_idx      ; cmat_ch_6(1,2) = 9.66e-3_idx   ; cmat_ch_6(1,3) = -2.28e-5_idx
    cmat_ch_6(2,1) = -2.67e-3_idx   ; cmat_ch_6(2,2) = 2.103e-1_idx  ; cmat_ch_6(2,3) = -5.329_idx
    cmat_ch_6(3,1) =  1.55e-3_idx   ; cmat_ch_6(3,2) = -6.228e-2_idx ; cmat_ch_6(3,3) =  5.094e-1_idx

    t_diff = t_a-t_s-(t_a+273.16_idx)*0.61_idx*(q_s-q_a)
  
    if (ws10 <= 5.0_idx) then
        if (t_diff <= -0.75_idx) then
           C_D0 = cmat_cd_1(1,1) + cmat_cd_1(1,2) * ws10 + cmat_cd_1(1,3) * (ws10**2) + cmat_cd_1(1,4) * (ws10**3)
           C_D1 = cmat_cd_1(2,1) + cmat_cd_1(2,2) / ws10 + cmat_cd_1(2,3) / (ws10**2) + cmat_cd_1(2,4) / (ws10**3)
           C_D2 = cmat_cd_1(3,1) + cmat_cd_1(3,2) / ws10 + cmat_cd_1(3,3) / (ws10**2) + cmat_cd_1(3,4) / (ws10**3)           
           C_L0 = cmat_ch_1(1,1) + cmat_ch_1(1,2) * ws10 + cmat_ch_1(1,3) * (ws10**2)
           C_L1 = cmat_ch_1(2,1) + cmat_ch_1(2,2) * ws10 + cmat_ch_1(2,3) * (ws10**2)
           C_L2 = cmat_ch_1(3,1) + cmat_ch_1(3,2) * ws10 + cmat_ch_1(3,3) * (ws10**2)
        else if (t_diff < 0.75_idx) then
           C_D0 = cmat_cd_2(1,1) + cmat_cd_2(1,2) * ws10 + cmat_cd_2(1,3) * (ws10**2)
           C_D1 = cmat_cd_2(2,1) + cmat_cd_2(2,2) / ws10 + cmat_cd_2(2,3) / (ws10**2)
           C_D2 = 0.0_idx
           C_L0 = cmat_ch_2(1,1) + cmat_ch_2(1,2) * ws10 + cmat_ch_2(1,3) * (ws10**2)
           C_L1 = cmat_ch_2(2,1) + cmat_ch_2(2,2) * ws10 + cmat_ch_2(2,3) * (ws10**2)
           C_L2 = 0.0_idx
        else
           C_L0 = cmat_ch_3(1,1) + cmat_ch_3(1,2) * ws10 + cmat_ch_3(1,3) * (ws10**2)
           C_L1 = cmat_ch_3(2,1) + cmat_ch_3(2,2) * ws10 + cmat_ch_3(2,3) * (ws10**2)
           C_L2 = cmat_ch_3(3,1) + cmat_ch_3(3,2) * ws10 + cmat_ch_3(3,3) * (ws10**2)           
        end if
     else
        if (t_diff <= -0.75_idx) then
           C_L0 = cmat_ch_4(1,1) + cmat_ch_4(1,2) * ws10 + cmat_ch_4(1,3) * (ws10**2)
           C_L1 = cmat_ch_4(2,1) + cmat_ch_4(2,2) / ws10 + cmat_ch_4(2,3) / (ws10**2)
           C_L2 = cmat_ch_4(3,1) + cmat_ch_4(3,2) / ws10 + cmat_ch_4(3,3) / (ws10**2)           
        else if (t_diff < 0.75) then
           C_L0 = cmat_ch_5(1,1) + cmat_ch_5(1,2) * ws10 + cmat_ch_5(1,3) * (ws10**2)
           C_L1 = cmat_ch_5(2,1) + cmat_ch_5(2,2) / ws10 + cmat_ch_5(2,3) / (ws10**2)
           C_L2 = 0.0_idx
        else
           C_L0 = cmat_ch_6(1,1) + cmat_ch_6(1,2) * ws10 + cmat_ch_6(1,3) * (ws10**2)
           C_L1 = cmat_ch_6(2,1) + cmat_ch_6(2,2) / ws10 + cmat_ch_6(2,3) / (ws10**2)
           C_L2 = cmat_ch_6(3,1) + cmat_ch_6(3,2) / ws10 + cmat_ch_6(3,3) / (ws10**2)           
        end if
     end if
     C_h = (C_L0 + C_L1 * t_diff + C_L2 * t_diff**2)*1.0e-3_idx
     C_e = C_h
     sh = -1.0_idx * rho_a * C_pa * C_h * ws10 * (t_s - t_a)
     lh = -1.0_idx * rho_a * L * C_e * ws10 * (q_s - q_a)
     ev = lh_to_ev * lh
!     ev = rho_a * C_e * ws10 * (q_s-q_a)*1.0e-3_idx
   end subroutine bulk_kara05
  subroutine bulk_ncep(t_a,q_a,t_s,ws,sh,lh,ev)
    implicit none
    real(idx),intent(in) :: t_a,q_a,t_s,ws
    real(idx),intent(inout) :: sh,lh,ev
    real(idx) :: t_a_K,t_s_K
    real(idx) :: p_air,R_g,C_pa,e_s,q_s,L,rho_a
    real(idx) :: u,u10,t10,q10
    real(idx) :: cd_n10,cd_n10_rt,ce_n10,stab,ch_n10
    real(idx) :: tv, ustar,bstar,tstar, qstar,z0,xx
    real(idx) :: cd,ch,ce,cd_rt
    integer, parameter :: n_itts = 10
    integer :: j
    real(idx), parameter :: GRAV   = 9.80
    real(idx), parameter :: VONKARM = 0.40
    real(idx) :: zetat,psi_mt, psi_ht            ! stability parameters
    real(idx) :: zetaq,psi_mq, psi_hq            ! stability parameters
    real(idx) :: zetaw, x2w, xw, psi_mw, psi_hw            ! stability parameters
    real(idx),parameter :: z_w=10.0,z_t = 2.0,z_q = 2.0 ! in [m]
    t_a_K=t_a+273.16_idx
    t_s_K=t_s+273.16_idx
    p_air = 1.013e3_idx ! hPa
    R_g = 2.871e2_idx  ! J kg^(-1) K^(-1)
    C_pa = 1.00467e3_idx  !  J kg^(-1) K (^-1)
    e_s=0.98_idx*((1.0007_idx+3.46e-6_idx*p_air)*6.1121_idx*exp(17.502_idx*t_s/(240.97+t_s)))
    q_s = 0.62197_idx * e_s / (p_air - 0.378_idx * e_s) ! [g/g]
    L = (1.0e6_idx) * (2.501 - 0.00237_idx * t_s) ! [J/kg]
    rho_a = (1.0e2_idx) * p_air / (R_g * (T_a + 273.16_idx) * (1.0_idx + 0.61_idx * q_a)) ! [kg/m^3]
    ! Bulk
    u = max(ws,0.5_idx);  ! 0.5 m/s floor on wind (undocumented NCAR)
    u10 = u; t10 = t_a_K ; q10 = q_a
    cd_n10 = (2.7_idx/u10+1.42e-1_idx+7.64e-2_idx*u10)/1.0e3_idx; ! L-Y eqn. 6a
    cd_n10_rt = sqrt(cd_n10);
    ce_n10 = 34.6_idx *cd_n10_rt/1.0e3_idx;       ! L-Y eqn. 6b
    stab = 0.5_idx + sign(0.5_idx,t_a-t_s)
    ch_n10 = (18.0_idx*stab+32.7_idx*(1.0_idx-stab))*cd_n10_rt/1.0e3_idx;! L-Y eqn. 6c
    cd = cd_n10; ch = ch_n10; ce = ce_n10
    do j=1,n_itts                                           ! Monin-Obukhov iteration
       tv = t10*(1+0.608*q10);
       cd_rt = sqrt(cd);
       ustar = cd_rt*u;                               ! L-Y eqn. 7a
       tstar    = (ch/cd_rt)*(t10-t_s_K);                ! L-Y eqn. 7b
       qstar    = (ce/cd_rt)*(q10-q_s);                ! L-Y eqn. 7c
       bstar = grav*(tstar/tv+qstar/(q_a+1/0.608));
       zetaw     = vonkarm*bstar*z_w/(ustar*ustar); ! L-Y eqn. 8a
       zetaw     = sign(min(abs(zetaw),10.0), zetaw);         ! undocumented NCAR
       zetat     = vonkarm*bstar*z_t/(ustar*ustar); ! L-Y eqn. 8a
       zetat     = sign(min(abs(zetat),10.0), zetat);         ! undocumented NCAR
       zetaq     = vonkarm*bstar*z_q/(ustar*ustar); ! L-Y eqn. 8a
       zetaq     = sign(min(abs(zetaq),10.0), zetaq);         ! undocumented NCAR

       call ncar_psi(zetaw,psi_mw,psi_hw)
       call ncar_psi(zetat,psi_mt,psi_ht)
       call ncar_psi(zetaq,psi_mq,psi_hq)
       ! Shift to 10m
       u10 = u/(1+cd_n10_rt*(log(z_w/10)-psi_mw)/vonkarm);       ! L-Y eqn. 9a
       t10 = t_a_K - tstar * (log(z_t/z_w)+psi_hw-psi_ht)/ vonkarm    ! L-Y eqn. 9
       q10 = q_a - qstar * (log(z_q/z_w)+psi_hw-psi_hq)/ vonkarm    ! L-Y eqn. 9
       cd_n10 = (2.7/u10+0.142+0.0764*u10)/1e3;                  ! L-Y eqn. 6a again
       cd_n10_rt = sqrt(cd_n10);
       ce_n10 = 34.6*cd_n10_rt/1e3;                              ! L-Y eqn. 6b again
       stab = 0.5 + sign(0.5_idx,zetaw)
       ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3;         ! L-Y eqn. 6c again
       z0 = 10*exp(-vonkarm/cd_n10_rt);                          ! diagnostic
       xx = (log(z_w/10)-psi_mw)/vonkarm;
       cd = cd_n10/(1+cd_n10_rt*xx)**2;                       ! L-Y 10a
       xx = (log(z_w/10)-psi_hw)/vonkarm;
       ch = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10) ! 10b (corrected code aug2007)
       ce = ce_n10/(1+ce_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10) ! 10c (corrected code aug2007)
    end do
    sh = -1.0_idx * rho_a * C_pa * ch * ws * (t_s - t_a)
    lh = -1.0_idx * rho_a * L * ce * ws * (q_s - q_a)
    ev = lh_to_ev * lh
!    ev = rho_a * ce * ws * (q_s-q_a) * 1.0e-3_idx
  end subroutine bulk_ncep
  subroutine ncar_psi(zeta,psim,psih)
     implicit none
     real(idx),intent(in) :: zeta
     real(idx) :: x,x2
     real(idx),intent(out) :: psim,psih
     x2 = sqrt(abs(1.0_idx-16.0_idx*zeta)) ! L-Y eqn. 8b
     x2 = max(x2, 1.0_idx);           ! undocumented NCAR
     x = sqrt(x2);
     if (zeta > 0) then
        psim = -5.0_idx*zeta; ! L-Y eqn. 8c
        psih = -5.0_idx*zeta; ! L-Y eqn. 8c
     else
        psim = log((1.0_idx+2.0_idx*x+x2)*(1.0_idx+x2)/8.0_idx)-2.0_idx*(atan(x)-atan(1.0_idx)); ! L-Y eqn. 8d
        psih = 2.0_idx*log((1.0_idx+x2)/2.0_idx);                                ! L-Y eqn. 8e
    end if
  end subroutine ncar_psi
   function diurnal_sw(SW,time) result(QSW)
    implicit none
    real(idx),intent(in) :: SW,time
    real(idx) :: QSW
    real(idx) :: time_24
    integer :: nday
    real(idx),parameter :: pi=4.0_idx * atan(1.0_idx)
    nday = time / (60.0_idx*60.0_idx*24.0_idx)
    time_24 =24.0_idx * (time / (60.0_idx*60.0_idx*24.0_idx) - nday)
    if (time_24 >= 6.0_idx .and. time_24 <= 18.0_idx) then
       QSW = pi * SW * sin(2.0_idx * pi * (time_24-6.0_idx)/24.0_idx)
    else
       QSW = 0.0_idx
    end if
  end function diurnal_sw

  function ld_to_lw(ld,sst) result(lw)
    implicit none
    real(idx),intent(in) :: ld,sst
    real(idx) :: lw
    real(idx) :: StefanBoltzman=5.67e-8
    lw = ld-StefanBoltzman*((sst+273.16_idx)**4)
  end function ld_to_lw
end module bulk_sub
