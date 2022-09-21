module nnf_sub
  use param
  use eos_sub
  implicit none
  private
  !-------------------------------------------------------  
  !Closure constants
  real(idx),parameter :: A_1=1.18_idx,A_2=0.665_idx
  real(idx),parameter :: B_1=24.0_idx,B_2=15.0_idx
  real(idx),parameter :: C_1 = 0.137_idx,C_2=0.75_idx,C_3=0.352_idx,C_5=0.2_idx
  !-------------------------------------------------------
  real(idx),parameter :: SQ=3.0_idx 
  real(idx),parameter :: q2_min = 1.0e-8_idx,q2l_min= 1.0e-8_idx
  real(idx),parameter :: gamma1=1.0_idx/3.0_idx-2.0_idx*A_1/B_1
  real(idx),parameter :: gamma2=(2.0_idx*A_1*(3.0_idx-2.0_idx*C_2)+B_2*(1.0_idx-C_3)) / B_1
  !F1
  real(idx),parameter :: F_1=B_1*(gamma1-C_1)+2.0_idx*A_1*(3.0_idx-2.0_idx*C_2)+&
       & 3.0_idx*A_2*(1.0_idx-C_2)*(1.0_idx-C_5)
  real(idx),parameter :: F_2=B_1*(gamma1+gamma2)-3.0_idx*A_1*(1.0_idx-C_2)
  !R_f1
  real(idx),parameter :: R_f1=B_1 * (gamma1-C_1) / F_1
  !R_f2
  real(idx),parameter :: R_f2=B_1*gamma1/F_2
  !R_fc: Critical flux Richardson number
  real(idx),parameter :: R_fc=gamma1 / (gamma1+gamma2)
  real(idx) :: tiny=1.0e-12_idx
  ! public setting--------------------------------------------------------
!  public :: cal_dens,cal_alpha,cal_beta,cal_bv,cal_shear2
!  public :: heat_absorb
!  public :: initialize_turb,initialize_turb2
  public :: cal_SH_my2,cal_SM_my2,cal_Rf
  public :: cal_SH_mynnf25,cal_SM_mynnf25
  public :: cal_mo_inv,diag_l_mynnf,mynnf25_vdiff,cal_kq
  public :: B_1
contains
  ! !============================================================================
  ! function cal_dens(s,t) result(ret)
  !   implicit none
  !   real(idx),intent(in) :: s,t
  !   real(idx) :: rho_w,ret
  !   rho_w = 999.842594_idx + t * 6.793952e-2_idx - (t**2) *  9.095290e-3_idx &
  !        &          + 1.001685e-04_idx * t**3 - 1.120083e-6_idx * t**4 + 6.536332e-9_idx * t**5 
  !   ret = rho_w  + (0.824493_idx-4.0899e-3_idx * t + (7.6438e-5_idx) * t**2 -8.2467e-7_idx &
  !        &               * t**3 + 5.3875e-9_idx * t**4)*s &
  !        &                + (-5.72466e-3_idx + 1.0227e-4_idx * t - 1.6546e-6_idx * t**2) * s**(1.5_idx) &
  !        &                + 4.8314e-4_idx * s**2
  ! end function cal_dens
  ! !============================================================================
  ! !-Calculate thermal expansion coefficient
  ! ! Note that this function generally returns negative value
  ! function cal_alpha(s,t) result(alpha)
  !   implicit none
  !   real(idx),intent(in) :: s,t
  !   real(idx) :: rho_1,rho_2,dt
  !   real(idx) :: alpha
  !   dt=0.01_idx
  !   rho_1=cal_dens(s,t)
  !   rho_2=cal_dens(s,t+dt)
  !   alpha = (rho_2-rho_1)/(rho_1*dt)  ! alpha < 0
  ! end function cal_alpha
  ! !============================================================================
  ! !-Calculate salinity extraction coefficient
  ! function cal_beta(s,t) result(beta)
  !   implicit none
  !   real(idx),intent(in) :: s,t
  !   real(idx) :: rho_1,rho_2,ds
  !   real(idx) :: beta
  !   ds=0.01_idx
  !   rho_1=cal_dens(s+ds,t)
  !   rho_2=cal_dens(s,t)    ! rho_1 > rho_2
  !   beta = (rho_1-rho_2)/(rho_1*ds)  ! beta > 0
  ! end function cal_beta
  ! !============================================================================
  ! function cal_bv(nz,z_rho,temp,salt) result(bv)
  !   implicit none
  !   integer,intent(in) :: nz
  !   real(idx),intent(in) :: z_rho(nz),temp(nz),salt(nz)
  !   real(idx) :: bv(nz-1)
  !   integer :: iz
  !   do iz=1,nz-1
  !      ! Buoyancy frequency
  !      bv(iz) =  - 1.0_idx * g * (cal_dens(salt(iz+1),temp(iz+1))- cal_dens(salt(iz),temp(iz))) / &
  !           & ((z_rho(iz+1)-z_rho(iz))* cal_dens(salt(iz),temp(iz)))
  !   end do
  ! end function cal_bv
  ! function cal_shear2(nz,z_rho,u,v) result(shear2)
  !   implicit none
  !   integer,intent(in) :: nz
  !   real(idx),intent(in) :: z_rho(nz),u(nz),v(nz)
  !   real(idx) :: shear2(nz-1)
  !   real(idx) :: uz,vz
  !   integer :: iz
  !   do iz=1,nz-1
  !      ! Shear
  !      uz = (u(iz+1)-u(iz)) / (z_rho(iz+1)-z_rho(iz))
  !      vz = (v(iz+1)-v(iz)) / (z_rho(iz+1)-z_rho(iz))
  !      shear2(iz) = (uz**2 + vz**2 + tiny)
  !   end do
  ! end function cal_shear2
  ! !================================================
  ! ! Initialize turbulence array
  ! subroutine initialize_turb(N,z,T,S,U,V,q2,l)
  !   implicit none
  !   integer,intent(in) :: N
  !   real(idx),intent(in) :: z(N),T(N),S(N),U(N),V(N)
  !   real(idx),intent(inout) :: q2(N),l(N)
  !   real(idx),parameter :: H=20.0_idx
  !   integer :: i
  !   ! initial condition=========================================
  !   do i=1,N
  !      q2(i) = q2_min!*exp(z(i)/H)
  !      l(i) = kappa * abs(z(i))
  !   end do
  ! end subroutine initialize_turb
  ! subroutine initialize_turb2(N,z,q2,l)
  !   implicit none
  !   integer,intent(in) :: N
  !   real(idx),intent(in) :: z(N)
  !   real(idx),intent(inout) :: q2(N),l(N)
  !   real(idx),parameter :: H=20.0_idx
  !   integer :: i
  !   ! initial condition=========================================
  !   do i=1,N
  !      q2(i) = q2_init*exp(z(i)/H)
  !      l(i) = kappa * abs(z(i))
  !   end do
  ! end subroutine initialize_turb2
  ! !============================================================================
  ! function cal_penet(level,sw,beta1,beta2,r_long) result(ret)
  !   implicit none
  !   real(idx) :: level,sw,beta1,beta2,r_long
  !   real(idx) :: ret
  !   real(idx) :: z1,z1b1,z1b2
  !   z1=abs(level)
  !   z1b1 = min(z1/beta1,100.0_idx) ; z1b2 = min(z1/beta2,100.0_idx)
  !   ret = sw * (r_long * exp(-1.0_idx * z1b1) +&
  !        & (1.0_idx - r_long) * exp(-1.0_idx * z1b2))
  ! end function cal_penet
  ! subroutine heat_absorb(nz,level,heat_sl,beta1,beta2,r_long,heat)
  !   implicit none
  !   integer,intent(in) :: nz
  !   real(idx),intent(in) :: level(nz)
  !   real(idx),intent(in) :: heat_sl
  !   real(idx),intent(in) :: beta1,beta2,r_long
  !   real(idx),intent(inout) :: heat(nz)
  !   integer :: i
  !   do i=1,nz
  !      heat(i) = cal_penet(level(i),heat_sl,beta1,beta2,r_long)
  !   end do
  ! end subroutine heat_absorb
  !================================================
  ! calculate flux richardson number 
  function cal_Rf(Ri) result(Rf)
    implicit none
    real(idx),intent(in) :: Ri
    real(idx) :: Rf
    real(idx) :: R_i1,R_i2,R_i3
    R_i1 = 0.5_idx * A_2 * F_2 / (A_1*F_1)
    R_i2 = 0.5_idx * R_f1 / R_i1
    R_i3 = (2.0_idx * R_f2 - R_f1) / R_i1
    Rf = R_i1 * (Ri + R_i2-sqrt(Ri**2-R_i3*Ri+R_i2**2))
    if (Rf .ge. R_fc) then
       Rf=R_fc-tiny
    end if
  end function cal_Rf
  function cal_SH_my2(Rf) result(SH2)
    implicit none
    real(idx),intent(in) :: Rf
    real(idx) :: SH2
    SH2 =3.0_idx  * A_2 * (gamma1+gamma2) * (R_fc-Rf) / (1.0_idx-Rf)
  end function cal_SH_my2
  function cal_SM_my2(Rf) result(SM2)
    implicit none
    real(idx),intent(in) :: Rf
    real(idx) :: SH2,SM2
    SH2 = cal_SH_my2(Rf)
    SM2 = ((A_1*F_1)/ (A_2 * F_2)) * ((R_f1-Rf)/(R_f2-Rf))*SH2
  end function cal_SM_my2
  !============================================================================
  ! Mellor-Yamada Nakanishi-Niino-Furuichi parameterization
  function cal_SH_mynnf25(GH,GM,alpha_C) result(SH)
    implicit none
    real(idx),intent(in) :: GH,GM,alpha_C
    real(idx) :: SH
    real(idx) :: denom,numer
    real(idx) :: phi_1,phi_2,phi_3,phi_4,phi_5
    phi_5 = 1.0_idx - 3.0_idx * (alpha_C**2) * A_2 * B_2 * (1.0_idx - C_3) * GH
    phi_1 = 1.0_idx - 9.0_idx * (alpha_C**2) * A_1 * A_2 * (1.0_idx - C_2) * GH
    phi_2 = phi_5 + 9.0_idx * (alpha_C**2) * (A_2**2) * (1.0_idx - C_2) * (1.0_idx-C_5) * GH
    phi_3 = phi_5 - 12.0_idx * (alpha_C**2) * (A_1*A_2) * (1.0_idx - C_2) * GH
    phi_4 = 6.0_idx * (alpha_C**2) * (A_1**2) * GM
    numer = alpha_C * A_2 * (phi_1 + 3.0_idx * C_1 * phi_4)
    denom = phi_1 * phi_3 + phi_2 * phi_4
    SH=numer/denom
  end function cal_SH_mynnf25
  function cal_SM_mynnf25(GH,GM,alpha_C) result(SM)
    implicit none
    real(idx),intent(in) :: GH,GM,alpha_C
    real(idx) :: SM
    real(idx) :: denom,numer
    real(idx) :: phi_1,phi_2,phi_3,phi_4,phi_5
    phi_5 = 1.0_idx - 3.0_idx * (alpha_C**2) * A_2 * B_2 * (1.0_idx - C_3) * GH
    phi_1 = 1.0_idx - 9.0_idx * (alpha_C**2) * A_1 * A_2 * (1.0_idx - C_2) * GH
    phi_2 = phi_5 + 9.0_idx * (alpha_C**2) * (A_2**2) * (1.0_idx - C_2) * (1.0_idx - C_5) * GH
    phi_3 = phi_5 - 12.0_idx * (alpha_C**2) * (A_1 * A_2) * (1.0_idx - C_2) * GH
    phi_4 = 6.0_idx * (alpha_C**2) * (A_1**2) * GM
    numer = alpha_C * A_1 * (phi_2 - 3.0_idx * C_1 * phi_3)
    denom = phi_1 * phi_3 + phi_2 * phi_4
    SM = numer / denom
  end function cal_SM_mynnf25
  !============================================================================
  function cal_mo_inv(sst,sss,tau_x,tau_y,sw,nsw,e_p,beta1,beta2,r_long) result(hmo_inv)
    real(idx),intent(in) :: sst,sss,tau_x,tau_y,sw,nsw,e_p,beta1,beta2,r_long
    real(idx) :: hmo_inv
    real(idx) :: alpha,beta,bf,tau,u_star2,u_star3
    real(idx) :: bf_nonsolar,bf_solar
    alpha = cal_alpha(sss,sst) ! < 0
    beta = cal_beta(sss,sst)   ! > 0
    !bf=-wb0=-g*alpha*wto+beta*ws0
    bf_nonsolar = -1.0_idx * g * (alpha * nsw / (rho * cpw)+ beta * S0 * e_p)
    bf_solar =-1.0_idx * g * alpha * sw / (rho * cpw)
    bf = bf_solar + bf_nonsolar
    ! tau_x = rho * u_star^2 * (direction)
    ! tau_y = rho * v_star^2
    ! u_star3 = (sqrt(u_star^2+v_star^2))^3
    u_star2=sqrt(tau_x**2+tau_y**2)/rho
    u_star3 = sqrt(u_star2)**3

    hmo_inv = kappa * bf / (max(u_star3,tiny))
    if (hmo_inv .gt. tiny) then
       bf_solar =-1.0_idx * g * alpha * (sw  - cal_penet(1.0_idx/max(hmo_inv,tiny),sw,beta1,beta2,r_long))/ (rho*cpw)

       bf = bf_solar + bf_nonsolar
       hmo_inv = kappa * bf / (max(u_star3,tiny))
    end if
    
    !bf <0 -> mo < 0 ! unstable
    !bf >0 -> mo > 0 ! stable
  end function cal_mo_inv
  !=========================================================
  ! Diagnoze length of eddy scale
  !=========================================================
  subroutine diag_l_mynnf(N,ZZ,bvf,Q2,L,tau_x,tau_y,sw,nsw,e_p,beta1,beta2,r_long)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: ZZ(N),bvf(N-1),Q2(N)
    real(idx),intent(in) :: tau_x,tau_y,sw,nsw,e_p,beta1,beta2,r_long
    real(idx),intent(inout) :: L(N)
    real(idx) :: mo_inv,xi
    real(idx) :: ls_inv,lt_inv,lb_inv,qz_int,q_int,l_inv
    integer :: iz
    mo_inv = cal_mo_inv(20.0_idx,32.0_idx,tau_x,tau_y,sw,nsw,e_p,beta1,beta2,r_long)
    L(N)=0.0_idx
    q_int=0.0_idx
    qz_int=0.0_idx
    do iz = n-1,1,-1
       ! Integration
       q_int=q_int + (zz(iz+1)-zz(iz))*(sqrt(Q2(iz+1))+sqrt(Q2(iz)))*0.5_idx
       qz_int=qz_int + (zz(iz+1)-zz(iz))*(sqrt(Q2(iz+1))*abs(zz(iz+1))+sqrt(Q2(iz))*abs(zz(iz)))*0.5_idx
    end do
    do iz=N-1,1,-1
       xi = abs(zz(iz)) * mo_inv
       ! Ls
       if (xi .ge. 1) then
          ls_inv=3.7_idx / (kappa*abs(zz(iz)))
       else if (xi .ge. 0) then
          ls_inv=(1.0_idx+2.7_idx*xi) / (kappa*abs(zz(iz)))
       else
          ls_inv=((1.0_idx-100.0_idx*xi)**(-0.2_idx)) / (kappa*abs(zz(iz)))
       end if
       ! Lt
       lt_inv= q_int / (0.23_idx*qz_int)
       ! Lb
       if (bvf(iz) .ge. 0) then
          lb_inv=sqrt(bvf(iz)/Q2(iz)) / 0.53_idx
       else
          lb_inv=0.0_idx+tiny
       end if
       l_inv = ls_inv+lt_inv+lb_inv
       !l_inv = lt_inv
       L(iz)=1.0_idx / l_inv
    end do
  end subroutine diag_l_mynnf
  !============================================================================
  ! NNFH subroutine
  subroutine mynnf25_vdiff(N,bvf,shear2,Q2,L,KM,KT,KS)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: bvf(N-1),shear2(N-1)
    real(idx),intent(inout) :: Q2(N),L(N),KM(N-1),KT(N-1),KS(N-1)
    real(idx) :: Ri,Rf,SM2,SH2,GH,GM,SM,SH,alpha_C
    integer :: i
    !-------------------------------------------------
    ! Calculate vertical differential　of each variable
    !-------------------------------------------------
    do i=1,N-1
       ! Level2 energy--------------------------
       ! GH (see p287 of Furuichi_et_al(2012))
       GH = - 1.0_idx * bvf(i) * l(i) * l(i) / q2(i)
       ! GM (see p287 of Furuichi_et_al(2012))
       GM = Shear2(i)  * l(i) * l(i) / q2(i)
       ! Richardson number
       Ri = bvf(i) / Shear2(i)
       ! Calculate Rf in (A11)
       Rf=cal_Rf(Ri)
       SM2=cal_SM_my2(Rf)
       alpha_C = sqrt(1.0_idx / max(tiny,(B_1*cal_SM_my2(Rf)*(1.0_idx-Rf)*GM)))
       alpha_C = min(alpha_C,1.0_idx)
      
       SM=cal_SM_mynnf25(GH,GM,alpha_C)
       SH=cal_SH_mynnf25(GH,GM,alpha_C)
       !------------------------------------------------
       ! Compute vertical viscosity
       !------------------------------------------------
       km(i) = l(i) * sqrt(q2(i)) * SM
       kt(i) = l(i) * sqrt(q2(i)) * SH
       ks(i) = l(i) * sqrt(q2(i)) * SH
    end do
  end subroutine mynnf25_vdiff
  ! calculate Kq
  subroutine cal_kq(n,km,kq)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: km(n-1)
    real(idx),intent(inout) :: kq(n-1)
    integer :: i
    do i =1,n-2
       kq(i) = 0.5_idx * (km(i)+km(i+1)) * Sq
    end do
    kq(n-1) = 0.5_idx * km(n-1) * Sq
  end subroutine cal_kq
  !============================================================================
end module nnf_sub

