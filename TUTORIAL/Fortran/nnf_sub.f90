module nnf_sub
  use ml_param
  use ml_utils
  implicit none
  private
  !-------------------------------------------------------  
  !Closure constants
  real(idx),parameter :: A_1=1.18_idx,A_2=0.665_idx
  real(idx),parameter :: B_1=24.0_idx,B_2=15.0_idx
  real(idx),parameter :: C_1 = 0.137_idx,C_2=0.75_idx,C_3=0.352_idx,C_5=0.2_idx
  !-------------------------------------------------------
  real(idx),parameter :: SQ=3.0_idx 
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
  public :: cal_SH_my2,cal_SM_my2,cal_Rf
  public :: cal_SH_mynnf25,cal_SM_mynnf25
  public :: diag_l_mynnf,mynnf25_vdiff,cal_kq
  public :: B_1
contains
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
  ! Mellor-Yamada Nakanishi-Niino-Furuichi parameteriation
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
  !=========================================================!
  ! Diagnoze length of eddy scale
  !=========================================================!
  subroutine diag_l_mynnf(N,ZZ,bvf,Q2,L,tau_x,tau_y,sw,nsw,salflux)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: ZZ(N+1),bvf(N-1),Q2(N+1)
    real(idx),intent(in) :: tau_x,tau_y,sw,nsw,salflux
    real(idx),intent(inout) :: L(N+1)
    real(idx) :: mo_inv,xi,bvf_targ
    real(idx) :: ls_inv,lt_inv,lb_inv,qz_int,q_int,l_inv
    integer :: i
    mo_inv = cal_mo_inv(T0,S0,tau_x,tau_y,sw,nsw,salflux)
    L(N+1)=0.0_idx
    q_int=0.0_idx
    qz_int=0.0_idx
    do i = 1,n
       ! Integration
       q_int=q_int + (zz(i+1)-zz(i))*(sqrt(Q2(i+1))+sqrt(Q2(i)))*0.5_idx
       qz_int=qz_int + (zz(i+1)-zz(i))*(sqrt(Q2(i+1))*abs(zz(i+1))+sqrt(Q2(i))*abs(zz(i)))*0.5_idx
    end do
    do i=1,N
       xi = abs(zz(i)) * mo_inv
       ! Ls
       if (xi .ge. 1) then
          ls_inv=3.7_idx / (kappa*abs(zz(i)))
       else if (xi .ge. 0) then
          ls_inv=(1.0_idx+2.7_idx*xi) / (kappa*abs(zz(i)))
       else
          ls_inv=((1.0_idx-100.0_idx*xi)**(-0.2_idx)) / (kappa*abs(zz(i)))
       end if
!       write(*,*) zz(i),xi,ls_inv
       ! Lt
       lt_inv= q_int / (0.23_idx*qz_int)
       ! Lb
       if (i==1) then
          bvf_targ=bvf(1)
       else
          bvf_targ=bvf(i-1)
       end if
       if (bvf_targ .ge. 0) then
          lb_inv=sqrt(bvf_targ/(Q2(i)+tiny)) / 0.53_idx
       else
          lb_inv=0.0_idx+tiny
       end if
       l_inv = ls_inv+lt_inv+lb_inv
       L(i)=1.0_idx / l_inv
    end do
  end subroutine diag_l_mynnf
  !============================================================================
  ! NNFH subroutine
  subroutine mynnf25_vdiff(N,bvf,shear2,Q2,L,KM,KT,KS)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: bvf(N-1),shear2(N-1)
    real(idx),intent(inout) :: Q2(N+1),L(N+1),KM(N-1),KT(N-1),KS(N-1)
    real(idx) :: Ri,Rf,SM2,GH,GM,SM,SH,alpha_C
    integer :: i
    !-------------------------------------------------
    ! Calculate vertical differential　of each variable
    !-------------------------------------------------
    do i=1,N-1
       ! Level2 energy--------------------------
       ! GH (see p287 of Furuichi_et_al(2012))
       GH = - 1.0_idx * bvf(i) * l(i+1) * l(i+1) / q2(i+1)
       ! GM (see p287 of Furuichi_et_al(2012))
       GM = Shear2(i)  * l(i+1) * l(i+1) / q2(i+1)
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
       km(i) = l(i+1) * sqrt(q2(i+1)) * SM
       kt(i) = l(i+1) * sqrt(q2(i+1)) * SH
       ks(i) = l(i+1) * sqrt(q2(i+1)) * SH
    end do
  end subroutine mynnf25_vdiff
  ! calculate Kq
  subroutine cal_kq(n,km,kq)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: km(n-1)
    real(idx),intent(inout) :: kq(n)
    integer :: i
    kq(1)=0.5*Sq*km(1)
    do i =1,n-2
       kq(i+1) = 0.5_idx * (km(i)+km(i+1)) * Sq
    end do
    kq(n) = km(n-1) * Sq
  end subroutine cal_kq
  !============================================================================
end module nnf_sub

