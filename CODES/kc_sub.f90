module kc_sub
  use eos_sub
  use param
  implicit none
  private
  !------------------------------------------
  real(idx),parameter :: A_1=0.92_idx,A_2=0.74_idx
  real(idx),parameter :: B_1=16.6_idx,B_2=10.1_idx
  real(idx),parameter :: C_1 = 0.08_idx,C_2=0.7_idx,C_3=0.2_idx,C_5=0.2_idx
  real(idx),parameter :: SQ=0.41_idx
  !------------------------------------------
  real(idx),parameter :: E_1=1.8_idx,E_2=1.33_idx,E_3=1.0_idx,E_4=1.0_idx
  real(idx),parameter :: q2_min = 1.0e-9_idx,q2l_min= 1.0e-9_idx
  real(idx),parameter :: q2_init = 1.0e-8_idx
  real(idx) :: tiny=1.0e-12_idx
  ! Public setting--------------------------------------------------------------
  public :: cal_kq
  public :: cal_SH_my25,cal_SM_my25,my25_vdiff
  public :: B_1,E_1,E_2
contains
  !============================================================================
  function cal_SM_my25(GH) result(SM)
    implicit none
    real(idx),intent(inout) :: GH
    real(idx) ::SH,SM
    real(idx) ::denom_1,denom_2,numer_1,numer_2
    if (GH .lt.-0.28_idx) then
       GH = -0.28_idx
    end if
    if (GH .ge.0.029_idx) then
       GH = 0.029_idx
    end if
    denom_1=1.0_idx-3.0_idx*A_2*GH*(6.0_idx*A_1+B_2*(1.0_idx-C_3))
    numer_1=A_2*(1.0_idx-6.0_idx*A_1/B_1)
    SH=numer_1/denom_1
    numer_2=A_1*((1.0_idx-6.0_idx*A_1/B_1-3.0_idx*C_1)+9.0_idx*(2.0_idx*A_1+A_2*(1.0_idx-C_2))*SH*GH)
    denom_2=1.0_idx-9.0_idx*A_1*A_2*GH
    SM=numer_2/denom_2
  end function cal_SM_my25
  function cal_SH_my25(GH) result(SH)
    implicit none
    real(idx),intent(in) :: GH
    real(idx) ::SH
    real(idx) ::denom,numer
    denom=1.0_idx-3.0_idx*A_2*GH*(6.0_idx*A_1+B_2*(1.0_idx-C_3))
    numer=A_2*(1.0_idx-6.0_idx*A_1/B_1)
    SH=numer/denom
  end function cal_SH_my25
  !============================================================================
  subroutine my25_vdiff(N,bvf,shear2,Q2,L,KM,KT,KS)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: bvf(N-1),shear2(N-1)
    real(idx),intent(inout) :: Q2(N),L(N),KM(N),KT(N),KS(N)
    real(idx) :: GH,SM,SH
    integer :: i
    !-------------------------------------------------
    ! Calculate vertical differential　of each variable
    !-------------------------------------------------
    do i=1,N-1
       !-------------------------------------------------
       ! Compute vertical viscocity coefficients
       !-------------------------------------------------
       ! GH (see Eq. (30) of KY94)
       GH= -1.0_idx * bvf(i) * l(i) * l(i) / q2(i)
       ! SH (see Eq. (28) of KY94)
       SH=cal_SH_my25(GH)
       ! SM (see Eq. (29) of KY94)
       SM=cal_SM_my25(GH)
       !------------------------------------------------
       ! Compute vertical viscosity
       !------------------------------------------------
       km(i) = l(i) * sqrt(q2(i)) * SM
       kt(i) = l(i) * sqrt(q2(i)) * SH
       ks(i) = l(i) * sqrt(q2(i)) * SH
    end do
  end subroutine my25_vdiff
  !============================================================================
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
    kq(n-1) =  0.5_idx * km(n-1) * Sq
  end subroutine cal_kq
end module kc_sub

