module solve_diag
  use param
  implicit none
  private
#include 'CPPLISTS.h'
  !-------------------------------------------------------
  real(idx),parameter :: qq_min = 1.0e-8_idx,qql_min= 1.0e-9_idx
  real(idx) :: tiny=1.0e-12_idx
  public :: solve_tri_implicit_real_nolap, solve_tri_implicit_complex_nolap
  public :: cal_next_theta,cal_next_theta_bud,cal_next_theta_cn_bud
  public :: cal_next_salt,cal_next_salt_bud,cal_next_salt_cn_bud
  public :: cal_next_uv,cal_next_uv_bud,cal_next_uv_cn_bud,cal_next_uv_dev_bud
  public :: cal_next_qq,cal_next_qq_bud,cal_next_qq_cn_bud
  public :: cal_next_l
contains
  !============================================================================
  subroutine solve_tri_implicit_real_nolap(N,A_in,B_in,C_in,D_in,U)
    ! Solve A * U_new(i-1) + B * U_new(i) + C * U_new(i+1) = D
    ! solve DT/dt=(k(z)(d^2 T)/(dz^2))
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: A_in(N),B_in(N),C_in(N),D_in(N)
    real(idx),intent(inout) :: U(N)
    real(idx) :: A(N),B(N),C(N),D(N)
    integer :: i
    real(idx) :: m
    A = A_in ; B = B_in ; C = C_in; D = D_in
    ! initialize array
    U(1:N) = 0.0_idx
    do i=2,N
       m = A(i) / B(i-1)
       B(i) = B(i) - m * C(i-1)
       D(i) = D(i) - m * D(i-1)
       U(i) = 0.0_idx
    end do
    U(N) = D(N) / B(N)
    do i = N-1,1,-1
       U(i) = (D(i)-C(i)*U(i+1)) / B(i)
    end do
  end subroutine solve_tri_implicit_real_nolap
  !============================================================================
  subroutine solve_tri_implicit_complex_nolap(N,A,B,C,D,U)
    ! Solve A * U_new(i-1) + B * U_new(i) + C * U_new(i+1) = D
    ! solve DT/dt=(k(z)(d^2 T)/(dz^2))
    implicit none
    integer,intent(in) :: N
    complex(kind(0d0)),intent(inout) :: A(N),B(N),C(N),D(N),U(N)
    integer :: i
    complex(kind(0d0)) :: m
    ! initialize array
    U(1:N) = 0.0_idx
    do i=2,N
       m = A(i) / B(i-1)
       B(i) = B(i) - m * C(i-1)
       D(i) = D(i) - m * D(i-1)
       U(i) = 0.0_idx
    end do
    U(N) = D(N) / B(N)
    do i = N-1,1,-1
       U(i) = (D(i)-C(i)*U(i+1)) / B(i)
    end do
  end subroutine solve_tri_implicit_complex_nolap

  !============================================================================
  ! subroutine for calculating new temperature
  subroutine cal_next_theta(dt,N,depth,zz,theta,kt,heat,heat_no_solar,heat_adv,theta_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),theta(N),kt(N-1),heat(N),heat_adv(N)
    real(idx),intent(inout) :: theta_next(N)
    real(idx),intent(in) :: heat_no_solar
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    integer :: i
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * dtdz_dw * Kt(i-1)      ! Coeff. of T(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * Kt(i)+dtdz_dw*Kt(i-1))    ! Coefficient of T(i)_new
       C(i) = - 1.0_idx * dtdz_up * Kt(i)                 ! Coefficient of T(i+1)_new
       D(i) = Theta(i) + dt * (heat(i)-heat(i-1)) / (rho*cpw*(zz(i)-zz(i-1)))  &
            & +dt * heat_adv(i)  ! Coefficient of RHS
    end do
    ! Boundary Condition----------------
    ! At the surface
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * dtdz_dw * Kt(N-1)       ! Coefficient of T(N-1)
    B(N) = 1.0_idx + dtdz_dw * Kt(N-1)        ! Coefficient of T(N)
    C(N) = 0.0_idx                         ! Coefficient of T(N+1)_new(does not exist)
    ! KT(N)*dT/dz = Q
    ! dT/dt = ((Q/rcph)-vdiff(N-1))+() /tau_rst
    D(N) =  Theta(N) + dt * heat_no_solar / (rho *cpw*(zz(N)-zz(N-1))) + &
         & dt * (heat(N)-heat(N-1)) / (rho*cpw*(zz(N)-zz(N-1)))+ dt * heat_adv(N)
    ! At the bottom
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(1) = 0.0_idx      ! Coeff. of T(0)_new
    B(1) = 1.0_idx + dtdz_up * Kt(1)    ! Coefficient of T(i)_new
    C(1) = - 1.0_idx * dtdz_up * Kt(1)                 ! Coefficient of T(i+1)_new
    D(1) = Theta(1)   ! Coefficient of RHS
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Theta_next)
  end subroutine cal_next_theta
  subroutine cal_next_theta_bud(dt,N,depth,zz,theta,kt,heat,heat_no_solar,heat_adv,theta_next,&
       & temp_tend,temp_vdiff,temp_penet,temp_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),theta(N),kt(N-1),heat(N),heat_adv(N)
    real(idx),intent(inout) :: theta_next(N)
    real(idx),intent(inout) :: temp_tend(N),temp_vdiff(N),temp_penet(N),temp_relax(N)
    real(idx),intent(in) :: heat_no_solar
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    integer :: i
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * dtdz_dw * Kt(i-1)      ! Coeff. of T(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * Kt(i)+dtdz_dw*Kt(i-1))    ! Coefficient of T(i)_new
       C(i) = - 1.0_idx * dtdz_up * Kt(i)                 ! Coefficient of T(i+1)_new
       D(i) = Theta(i) + dt * (heat(i)-heat(i-1)) / (rho*cpw*(zz(i)-zz(i-1)))  &
            & +dt * heat_adv(i)  ! Coefficient of RHS
    end do
    ! Boundary Condition----------------
    ! At the surface
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * dtdz_dw * Kt(N-1)       ! Coefficient of T(N-1)
    B(N) = 1.0_idx + dtdz_dw * Kt(N-1)        ! Coefficient of T(N)
    C(N) = 0.0_idx                         ! Coefficient of T(N+1)_new(does not exist)
    ! KT(N)*dT/dz = Q
    ! dT/dt = ((Q/rcph)-vdiff(N-1))+() /tau_rst
    D(N) =  Theta(N) + dt * heat_no_solar / (rho *cpw*(zz(N)-zz(N-1))) + &
         & dt * (heat(N)-heat(N-1)) / (rho*cpw*(zz(N)-zz(N-1)))+ dt * heat_adv(N)
    ! At the bottom
#ifdef zero_bottom
    A(1) = 0.0_idx      ! Coeff. of T(0)_new
    B(1) = 1.0_idx      ! Coefficient of T(1)_new
    C(1) = - 1.0_idx    ! Coefficient of T(1)_new
    D(1) = 0.0_idx      ! Ensure that T(1)=T(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(1) = 0.0_idx      ! Coeff. of T(0)_new
    B(1) = 1.0_idx + dtdz_up * Kt(1)    ! Coefficient of T(i)_new
    C(1) = - 1.0_idx * dtdz_up * Kt(1)                 ! Coefficient of T(i+1)_new
    D(1) = theta(1)+ dt * heat_adv(1)   ! Coefficient of RHS
#endif
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Theta_next)
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))

       temp_tend(i) = temp_tend(i) + (Theta_next(i) - Theta(i)) / dt
       temp_vdiff(i) = temp_vdiff(i) + Kt(i) * dtdz_up  * (Theta_next(i+1)-Theta_next(i)) - &
            & Kt(i-1) * dtdz_dw  * (Theta_next(i)-Theta_next(i-1))
       temp_penet(i) = temp_penet(i) + (heat(i)-heat(i-1)) / (rho*cpw*(zz(i)-zz(i-1)))
       temp_relax(i) = temp_relax(i) + heat_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    temp_tend(i) =  temp_tend(i)+ (Theta_next(i) - Theta(i)) / dt
    temp_vdiff(i) = temp_vdiff(i) + heat_no_solar / (rho *cpw*(zz(N)-zz(N-1))) &
         & - Kt(i-1) * dtdz_dw  * (Theta_next(i)-Theta_next(i-1))
    temp_penet(i) = temp_penet(i)+ (heat(N)-heat(N-1)) / (rho*cpw*(zz(N)-zz(N-1))) 
    temp_relax(i) = temp_relax(i) + heat_adv(i)
    i=1
#ifdef zero_bottom
    temp_tend(i) = temp_tend(i+1); temp_vdiff(i) = temp_vdiff(i+1)
    temp_penet(i) = temp_penet(i+1);  temp_relax(i) = temp_relax(i+1)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    temp_tend(i) =  temp_tend(i)+ (Theta_next(i) - theta(i)) / dt
    temp_vdiff(i) = temp_vdiff(i) + Kt(i) * dtdz_up  * (theta_next(i+1)-theta_next(i))
    temp_penet(i) = temp_penet(i)+ 0.0_idx
    temp_relax(i) = temp_relax(i) + heat_adv(i)
#endif
    !==================================
  end subroutine cal_next_theta_bud  

  subroutine cal_next_theta_cn_bud(dt,N,depth,zz,theta,kt,heat,heat_no_solar,heat_adv,theta_next,&
       & temp_tend,temp_vdiff,temp_penet,temp_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),theta(N),kt(N-1),heat(N),heat_adv(N)
    real(idx),intent(inout) :: theta_next(N)
    real(idx),intent(inout) :: temp_tend(N),temp_vdiff(N),temp_penet(N),temp_relax(N)
    real(idx),intent(in) :: heat_no_solar
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    integer :: i
    real(idx),parameter :: p=0.5_idx,q=1.0_idx-p
    real(idx) :: zz0
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * p * dtdz_dw * Kt(i-1)      ! Coeff. of T(i-1)_new
       B(i) = 1.0_idx + p *(dtdz_up * Kt(i)+dtdz_dw*Kt(i-1))    ! Coefficient of T(i)_new
       C(i) = - 1.0_idx * p * dtdz_up * Kt(i)                 ! Coefficient of T(i+1)_new
       D(i) = theta(i) + q * kt(i-1) * dtdz_dw * (theta(i-1)-theta(i)) &
            & + q * kt(i) * dtdz_up * (theta(i+1)-theta(i)) &
            & + dt * heat_adv(i)  &
            & + dt * (heat(i)-heat(i-1)) / (rho*cpw*(zz(i)-zz(i-1)))  ! Coefficient of RHS
    end do
    ! Boundary Condition----------------
    ! At the surface
    i=N
    dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    A(i) = - 1.0_idx * p * dtdz_dw * kt(i-1)      ! Coeff. of T(i-1)_new
    B(i) = 1.0_idx + p*dtdz_dw*kt(i-1)    ! Coefficient of T(i)_new
    C(i) = 0.0_idx
    D(i) = theta(i) + q * kt(i-1) * dtdz_dw * (theta(i-1)- theta(i)) &
         & + dt * heat_adv(i) &
         & + dt * (heat(i)-heat(i-1)) / (rho*cpw*(zz(i)-zz(i-1)))  &
         & + dt * heat_no_solar / (rho *cpw*(zz(N)-zz(N-1))) ! Coefficient of RHS
    ! At the bottom
    i=1
#ifdef zero_bottom
    A(i) = 0.0_idx      ! Coeff. of T(0)_new
    B(i) = 1.0_idx      ! Coefficient of T(1)_new
    C(i) = - 1.0_idx    ! Coefficient of T(1)_new
    D(i) = 0.0_idx      ! Ensure that T(1)=T(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(i) = 0.0_idx
    B(i) = 1.0_idx + p * dtdz_up * Kt(i)    ! Coefficient of T(i)_new
    C(i) = - 1.0_idx * p * dtdz_up * Kt(i)                 ! Coefficient of T(i+1)_new
    D(i) = theta(i) + q * kt(i) * dtdz_up * (theta(i+1)- theta(i)) &
         & + dt * heat_adv(i)  ! Coefficient of RHS
#endif
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Theta_next)
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       temp_tend(i) = temp_tend(i) + (Theta_next(i) - Theta(i)) / dt
       temp_vdiff(i) = temp_vdiff(i) + &
            & kt(i) * dtdz_up * (p*theta_next(i+1)+q*theta(i+1)-p*theta_next(i)-q*theta(i)) - &
            & kt(i-1) * dtdz_dw * (p*theta_next(i)+q*theta(i)-p*theta_next(i-1)-q*theta(i-1))
       temp_penet(i) = temp_penet(i) + (heat(i)-heat(i-1)) / (rho*cpw*(zz(i)-zz(i-1)))
       temp_relax(i) = temp_relax(i) + heat_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    temp_tend(i) =  temp_tend(i)+ (Theta_next(i) - Theta(i)) / dt
    temp_vdiff(i) = temp_vdiff(i) + heat_no_solar / (rho *cpw*(zz(N)-zz(N-1))) &
         & - kt(i-1) * dtdz_dw  * (p*theta_next(i)+q*theta(i)-p*theta_next(i-1)-q*theta_next(i-1))
    temp_penet(i) = temp_penet(i)+ (heat(N)-heat(N-1)) / (rho*cpw*(zz(N)-zz(N-1))) 
    temp_relax(i) = temp_relax(i) + heat_adv(i)
    i=1
#ifdef zero_bottom
    temp_tend(i) = temp_tend(i+1); temp_vdiff(i) = temp_vdiff(i+1)
    temp_penet(i) = temp_penet(i+1);  temp_relax(i) = temp_relax(i+1)
#else
    dtdz_up = 1.0_idx / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    temp_tend(i) =  temp_tend(i)+ (theta_next(i) - theta(i)) / dt
    temp_vdiff(i) = temp_vdiff(i) + kt(i) * dtdz_up * (p*theta_next(i+1)+q*theta(i+1)-p*theta_next(i)-q*theta(i))
    temp_penet(i) = temp_penet(i)+ 0.0_idx
    temp_relax(i) = temp_relax(i) + heat_adv(i)
#endif
    !==================================
  end subroutine cal_next_theta_cn_bud

  !---------------------------------
  ! Salinity
  !---------------------------------
  ! subroutine for calculating new salinity
  subroutine cal_next_salt(dt,N,depth,zz,salt,ks,ssflux,salt_adv,salt_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),salt(N),ks(N-1),salt_adv(N),ssflux
    real(idx),intent(inout) :: salt_next(N)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    integer :: i
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix-------------------------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * dtdz_dw * Ks(i-1)               ! Coefficient of S(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * Ks(i)+ dtdz_dw * Ks(i-1))         ! Coefficient of S(i)_new
       C(i) = - 1.0_idx * dtdz_up * Ks(i)                 ! Coefficient of S(i+1)_new
       D(i) = Salt(i) + dt * salt_adv(i)         ! Coefficient of RHS
    end do
    ! Boundary Condition--------------------------------
    ! At the surface---------------
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * dtdz_dw * Ks(N-1)    ! Coefficient of S(N-1)
    B(N) = 1.0_idx + dtdz_dw * Ks(N-1)     ! Coefficient of S(N)
    C(N) = 0.0_idx
    D(N) = Salt(N) + dt * ssflux / (zz(N)-zz(N-1)) + dt * salt_adv(N) ! KS(N-1)*dS/dz = S(N) * (E-p) => KS(N-1)*(S(N)-S(N-1)) / dz = (E-P)*S(N)
    ! At the bottom----------------
    A(1) = 0.0_idx     ! Coefficient of S(0)
    B(1) = -1.0_idx*Ks(1)/ (zz(2)-zz(1))          ! Coefficient of S(1)
    C(1) = Ks(1) / (zz(2)-zz(1))
    D(1) =  0.0_idx            ! KS(1)*dS/dz=0.0
    !---------------------------------------------------
    ! Solve
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Salt_next)
  end subroutine cal_next_salt
  subroutine cal_next_salt_bud(dt,N,depth,zz,salt,ks,ssflux,salt_adv,salt_next,&
       & salt_tend,salt_vdiff,salt_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),salt(N),ks(N-1),salt_adv(N),ssflux
    real(idx),intent(inout) :: salt_tend(N),salt_vdiff(N),salt_relax(N)
    real(idx),intent(inout) :: salt_next(N)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    integer :: i
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix-------------------------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * dtdz_dw * Ks(i-1)               ! Coefficient of S(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * Ks(i)+ dtdz_dw * Ks(i-1))         ! Coefficient of S(i)_new
       C(i) = - 1.0_idx * dtdz_up * Ks(i)                 ! Coefficient of S(i+1)_new
       D(i) = Salt(i) + dt * salt_adv(i)         ! Coefficient of RHS
    end do
    ! Boundary Condition--------------------------------
    ! At the surface---------------
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * dtdz_dw * Ks(N-1)    ! Coefficient of S(N-1)
    B(N) = 1.0_idx + dtdz_dw * Ks(N-1)     ! Coefficient of S(N)
    C(N) = 0.0_idx
    D(N) = salt(N) + dt * ssflux / (zz(N)-zz(N-1)) + dt * salt_adv(N) ! KS(N-1)*dS/dz = S(N) * (E-p) => KS(N-1)*(S(N)-S(N-1)) / dz = (E-P)*S(N)
        ! At the bottom----------------
#ifdef zero_bottom
    A(1) = 0.0_idx      ! Coeff. of S(0)_new
    B(1) = 1.0_idx      ! Coefficient of S(1)_new
    C(1) = - 1.0_idx    ! Coefficient of S(1)_new
    D(1) = 0.0_idx      ! Ensure that S(1)=S(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(1) = 0.0_idx      ! Coeff. of S(0)_new
    B(1) = 1.0_idx + dtdz_up * Ks(1)    ! Coefficient of S(1)_new
    C(1) = - 1.0_idx * dtdz_up * Ks(1)  ! Coefficient of S(2)_new
    D(1) = salt(1)+ dt * salt_adv(1)   ! Coefficient of RHS
#endif
    !---------------------------------------------------
    ! Solve
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Salt_next)
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))

       salt_tend(i) = salt_tend(i) + (Salt_next(i) - Salt(i)) / dt
       salt_vdiff(i) = salt_vdiff(i) + Ks(i) * dtdz_up  * (Salt_next(i+1)-Salt_next(i)) - &
            & Ks(i-1) * dtdz_dw  * (Salt_next(i)-Salt_next(i-1))
       salt_relax(i) = salt_relax(i) + salt_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    salt_tend(i) = salt_tend(i) + (Salt_next(i) - Salt(i)) / dt
    salt_vdiff(i) = salt_vdiff(i) + ssflux  / (zz(N)-zz(N-1)) - &
         & Ks(i-1) * dtdz_dw  * (Salt_next(i)-Salt_next(i-1))
    salt_relax(i) = salt_relax(i) + salt_adv(i)
    i=1
#ifdef zero_bottom
    salt_tend(i) = salt_tend(i+1); salt_vdiff(i) = salt_vdiff(i+1)
    salt_relax(i) = salt_relax(i+1)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    salt_tend(i) =  salt_tend(i)+ (salt_next(i) - salt(i)) / dt
    salt_vdiff(i) = salt_vdiff(i) + ks(i) * dtdz_up  * (salt_next(i+1)-salt_next(i))
    salt_relax(i) = salt_relax(i) + salt_adv(i)
#endif
  end subroutine cal_next_salt_bud

  subroutine cal_next_salt_cn_bud(dt,N,depth,zz,salt,ks,ssflux,salt_adv,salt_next,&
       & salt_tend,salt_vdiff,salt_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),salt(N),ks(N-1),salt_adv(N),ssflux
    real(idx),intent(inout) :: salt_tend(N),salt_vdiff(N),salt_relax(N)
    real(idx),intent(inout) :: salt_next(N)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    integer :: i
    real(idx),parameter :: p=0.5_idx,q=1.0_idx-p
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix-------------------------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * p * dtdz_dw * ks(i-1)               ! Coefficient of S(i-1)_new
       B(i) = 1.0_idx + p * (dtdz_up * ks(i)+ dtdz_dw * ks(i-1))         ! Coefficient of S(i)_new
       C(i) = - 1.0_idx * p * dtdz_up * ks(i)                 ! Coefficient of S(i+1)_new
       D(i) = salt(i) + q * dtdz_dw * ks(i-1)*(salt(i-1)-salt(i)) &
            & + q * dtdz_up * ks(i) * (salt(i+1)-salt(i)) &
            & + dt * salt_adv(i)        ! Coefficient of RHS
    end do
    ! Boundary Condition--------------------------------
    ! At the surface---------------
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * p * dtdz_dw * ks(N-1)    ! Coefficient of S(N-1)
    B(N) = 1.0_idx + p * dtdz_dw * ks(N-1)     ! Coefficient of S(N)
    C(N) = 0.0_idx
    D(N) = salt(N) + q * dtdz_dw * ks(N-1) * (salt(N-1) - salt(N)) &
         & + dt * ssflux / (zz(N)-zz(N-1)) + dt * salt_adv(N) ! KS(N-1)*dS/dz = S(N) * (E-p) => KS(N-1)*(S(N)-S(N-1)) / dz = (E-P)*S(N)
        ! At the bottom----------------
#ifdef zero_bottom
    A(1) = 0.0_idx      ! Coeff. of S(0)_new
    B(1) = 1.0_idx      ! Coefficient of S(1)_new
    C(1) = - 1.0_idx    ! Coefficient of S(1)_new
    D(1) = 0.0_idx      ! Ensure that S(1)=S(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(1) = 0.0_idx      ! Coeff. of S(0)_new
    B(1) = 1.0_idx + p * dtdz_up * ks(1)    ! Coefficient of S(1)_new
    C(1) = - 1.0_idx * p* dtdz_up * ks(1)  ! Coefficient of S(2)_new
    D(1) = salt(1) + q * dtdz_up * ks(1) * (salt(2)-salt(1)) &
         & + dt * salt_adv(1)   ! Coefficient of RHS
#endif
    !---------------------------------------------------
    ! Solve
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Salt_next)
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))

       salt_tend(i) = salt_tend(i) + (salt_next(i) - Salt(i)) / dt
       salt_vdiff(i) = salt_vdiff(i) + ks(i) * dtdz_up &
            & * (p*salt_next(i+1)+q*salt(i+1)-p*salt_next(i)-q*salt(i)) &
            & - ks(i-1) * dtdz_dw *(p*salt_next(i)+q*salt(i)-p*salt_next(i-1)-q*salt(i-1))
       salt_relax(i) = salt_relax(i) + salt_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    salt_tend(i) = salt_tend(i) + (Salt_next(i) - Salt(i)) / dt
    salt_vdiff(i) = salt_vdiff(i) + ssflux  / (zz(N)-zz(N-1)) - &
         & ks(i-1) * dtdz_dw  * (p*salt_next(i)+q*salt(i)-p*salt_next(i-1)-q*salt(i-1))
    salt_relax(i) = salt_relax(i) + salt_adv(i)
    i=1
#ifdef zero_bottom
    salt_tend(i) = salt_tend(i+1); salt_vdiff(i) = salt_vdiff(i+1)
    salt_relax(i) = salt_relax(i+1)
#else
    dtdz_up = 1.0_idx / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    salt_tend(i) =  salt_tend(i)+ (salt_next(i) - salt(i)) / dt
    salt_vdiff(i) = salt_vdiff(i) + ks(i) * dtdz_up  * (p*salt_next(i+1)+q*salt(i)-p*salt_next(i)-q*salt(i))
    salt_relax(i) = salt_relax(i) + salt_adv(i)
#endif
  end subroutine cal_next_salt_cn_bud

  !=================================================================
  !---------------------------------
  ! Velocity
  !---------------------------------
  ! subroutine for calculating new velocity
  subroutine cal_next_uv(dt,N,depth,zz,u,v,km,f,tau_x,tau_y,u_adv,v_adv,u_next,v_next)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: dt,depth(n),zz(n),u(n),v(n),km(n-1)
    real(idx),intent(in) :: f,tau_x,tau_y
    real(idx),intent(in) :: u_adv(n),v_adv(n)
    real(idx),intent(inout) :: u_next(n),v_next(n)
    complex(kind(0d0)) :: uv(n),uv_adv(n),A(N),B(N),C(N),D(N),uv_next(N)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i
    uv=0.0_idx ; uv_next = 0.0_idx; uv_adv=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       uv(i) = u(i) + (0.0_idx,1.0_idx) * v(i)
       uv_adv(i) = u_adv(i) + (0.0_idx,1.0_idx) * v_adv(i)
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = cmplx(-1.0_idx * dtdz_dw * Km(i-1),0.0_idx)               ! Coefficient of U(i-1)_new
       B(i) = cmplx(1.0_idx + (dtdz_dw*Km(i-1)+dtdz_up*Km(i)),f*dt)        ! Coefficient of U(i)_new
       C(i) = cmplx(-1.0_idx * dtdz_up * Km(i),0.0_idx)             ! Coefficient of U(i+1)_new
       D(i) = uv(i)+ dt * uv_adv(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    uv(N) = u(N) + (0.0_idx,1.0_idx) * v(N)
    uv_adv(N) = u_adv(N) + (0.0_idx,1.0_idx) * v_adv(N)

    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = cmplx(-1.0_idx * dtdz_dw * Km(N-1),0.0_idx)   ! Coefficient of U,V(N-1)        
    B(N) = cmplx(1.0_idx + dtdz_dw * Km(N-1),f*dt)       ! Coefficient of U,V(N)
    C(N) = cmplx(0.0_idx,0_idx)                   ! Coefficient of U,V(N+1)<-does not exist
    D(N) = uv(N) +  dt * cmplx(tau_x/rho,tau_y/rho) / (zz(N)-zz(N-1))+ dt * uv_adv(N)
    ! At the bottom
    A(1) = cmplx(0.0_idx,0.0_idx)                              ! Coefficient of U,V(0)<-does not exist
    B(1) = cmplx(-1.0_idx*Km(1) / (zz(2)-zz(1)),0.0_idx)       ! Coefficient of U,V(1)
    C(1) = cmplx(Km(1) / (zz(2)-zz(1)),0.0_idx)            ! Coefficient of U,V(2)
    D(1) = cmplx(0.0_idx,0.0_idx)
    call solve_tri_implicit_complex_nolap(N,A,B,C,D,uv_next)
    u_next=real(uv_next) ; v_next=imag(uv_next)

  end subroutine cal_next_uv
  subroutine cal_next_uv_bud(dt,N,depth,zz,u,v,km,f,tau_x,tau_y,u_adv,v_adv,u_next,v_next,&
       & u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),u(N),v(N),km(N-1)
    real(idx),intent(in) :: f,tau_x,tau_y
    real(idx),intent(in) :: u_adv(N),v_adv(N)
    real(idx),intent(inout) :: u_next(N),v_next(N)
    real(idx),intent(inout) :: u_rate(n),u_cor(n),u_vdiff(n),u_relax(n)
    real(idx),intent(inout) :: v_rate(n),v_cor(n),v_vdiff(n),v_relax(n)

    complex(kind(0d0)) :: uv(N),uv_adv(N),A(N),B(N),C(N),D(N),uv_next(N)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i
    uv=0.0_idx ; uv_next = 0.0_idx; uv_adv=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       uv(i) = u(i) + (0.0_idx,1.0_idx) * v(i)
       uv_adv(i) = u_adv(i) + (0.0_idx,1.0_idx) * v_adv(i)
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = cmplx(-1.0_idx * dtdz_dw * Km(i-1),0.0_idx)               ! Coefficient of U(i-1)_new
       B(i) = cmplx(1.0_idx + (dtdz_dw*Km(i-1)+dtdz_up*Km(i)),f*dt)        ! Coefficient of U(i)_new
       C(i) = cmplx(-1.0_idx * dtdz_up * Km(i),0.0_idx)             ! Coefficient of U(i+1)_new
       D(i) = uv(i)+ dt * uv_adv(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    uv(N) = u(N) + (0.0_idx,1.0_idx) * v(N)
    uv_adv(N) = u_adv(N) + (0.0_idx,1.0_idx) * v_adv(N)

    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(N) = cmplx(-1.0_idx * dtdz_dw * Km(N-1),0.0_idx)   ! Coefficient of U,V(N-1)        
    B(N) = cmplx(1.0_idx + dtdz_dw * Km(N-1),f*dt)       ! Coefficient of U,V(N)
    C(N) = cmplx(0.0_idx,0.0_idx)                   ! Coefficient of U,V(N+1)<-does not exist
    D(N) = uv(N) +  dt * cmplx(tau_x/rho,tau_y/rho) / (zz(N)-zz(N-1))+ dt * uv_adv(N)
    ! At the bottom
    i=1
#ifdef zero_bottom
    A(1) = cmplx(0.0_idx,0.0_idx)  ! Coefficient of U(0),V(0)<-does not exist
    B(1) = cmplx(-1.0_idx,0.0_idx) ! Coefficient of U(1),V(1)
    C(1) = cmplx(1.0_idx,0.0_idx)  ! Coefficient of U(2),V(2)
    D(1) = cmplx(0.0_idx,0.0_idx)  ! Ensure that U(1)=U(2); V(1)=V(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(1) = cmplx(0.0_idx,0.0_idx)  ! Coefficient of U(0),V(0)<-does not exist
    B(1) = cmplx(1.0_idx + dtdz_up*km(i),f*dt)        ! Coefficient of U(i)_new
    C(1) = cmplx(-1.0_idx * dtdz_up*km(i),0.0_idx)             ! Coefficient of U(i+1)_new
    D(1) = uv(1)+ dt * uv_adv(1)               ! Coefficient of RHS
#endif
    call solve_tri_implicit_complex_nolap(N,A,B,C,D,uv_next)
    u_next=real(uv_next) ; v_next=imag(uv_next)
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))       
       u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
       v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
       u_cor(i) = u_cor(i) + f * v_next(i)
       v_cor(i) = v_cor(i) - f * u_next(i)
       u_vdiff(i) = u_vdiff(i) + Km(i) * dtdz_up  * (u_next(i+1)-u_next(i)) - &
            & Km(i-1) * dtdz_dw  * (u_next(i)-u_next(i-1))
       v_vdiff(i) = v_vdiff(i) + Km(i) * dtdz_up  * (v_next(i+1)-v_next(i)) - &
            & Km(i-1) * dtdz_dw  * (v_next(i)-v_next(i-1))
       u_relax(i) = u_relax(i) + u_adv(i)
       v_relax(i) = v_relax(i) + v_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
    v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
    u_cor(i) = u_cor(i) + f * v_next(i)
    v_cor(i) = v_cor(i) - f * u_next(i)
    u_vdiff(i) = u_vdiff(i) + (tau_x/rho) / (zz(N)-zz(N-1)) - &
         & Km(i-1) * dtdz_dw  * (u_next(i)-u_next(i-1))
    v_vdiff(i) = v_vdiff(i) + (tau_y/rho) / (zz(N)-zz(N-1)) - &
         & Km(i-1) * dtdz_dw  * (v_next(i)-v_next(i-1))
    u_relax(i) = u_relax(i) + u_adv(i)
    v_relax(i) = v_relax(i) + v_adv(i)
    i=1
#ifdef zero_bottom
    u_rate(i) = u_rate(i+1); u_cor(i) = u_cor(i+1)
    u_vdiff(i) = u_vdiff(i+1); u_relax(i) = u_relax(i+1)
    v_rate(i) = v_rate(i+1);  v_cor(i) = v_cor(i+1)
    v_vdiff(i) = v_vdiff(i+1); v_relax(i) = v_relax(i+1)
#else
    dtdz_up = 1.0_idx / ((zz(2)-zz(1))*(depth(i+1)-depth(i)))       
    u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
    v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
    u_cor(i) = u_cor(i) + f * v_next(i); v_cor(i) = v_cor(i) - f * u_next(i)
    u_vdiff(i) = u_vdiff(i) + km(i) * dtdz_up  * (u_next(i+1)-u_next(i))
    v_vdiff(i) = v_vdiff(i) + km(i) * dtdz_up  * (v_next(i+1)-v_next(i))
    u_relax(i) = u_relax(i) + u_adv(i); v_relax(i) = v_relax(i) + v_adv(i)
#endif
  end subroutine cal_next_uv_bud

  subroutine cal_next_uv_cn_bud(dt,N,depth,zz,u,v,km,f,tau_x,tau_y,u_adv,v_adv,u_next,v_next,&
       & u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),u(N),v(N),km(N-1)
    real(idx),intent(in) :: f,tau_x,tau_y
    real(idx),intent(in) :: u_adv(N),v_adv(N)
    real(idx),intent(inout) :: u_next(N),v_next(N)
    real(idx),intent(inout) :: u_rate(n),u_cor(n),u_vdiff(n),u_relax(n)
    real(idx),intent(inout) :: v_rate(n),v_cor(n),v_vdiff(n),v_relax(n)
    complex(kind(0d0)) :: uv(N),uv_adv(N),A(N),B(N),C(N),D(N),uv_next(N)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i
    real(idx),parameter :: p=0.5_idx,q=1.0_idx-p,p1=0.0_idx,q1=0.0_idx-p1
    uv=0.0_idx ; uv_next = 0.0_idx; uv_adv=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       uv(i) = u(i) + cmplx(0.0_idx,1.0_idx) * v(i)
       uv_adv(i) = u_adv(i) + cmplx(0.0_idx,1.0_idx) * v_adv(i)
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = cmplx(-1.0_idx * p*dtdz_dw * km(i-1),0.0_idx)               ! Coefficient of U(i-1)_new
       B(i) = cmplx(1.0_idx + p*dtdz_dw*km(i-1)+p*dtdz_up*km(i),p1*f*dt)  ! Coefficient of U(i)_new
       C(i) = cmplx(-1.0_idx * p* dtdz_up * km(i),0.0_idx)             ! Coefficient of U(i+1)_new
       D(i) = uv(i) + q * dtdz_dw * km(i-1)*(uv(i-1)-uv(i)) &
            & + q * dtdz_up * km(i)*(uv(i+1)-uv(i)) !- dt * cmplx(0.0_idx,f) *q1 * uv(i)!&
       D(i) = uv(i) * (1.0_idx - q * dtdz_up * km(i)-q * dtdz_dw * km(i-1)) &
            + q * dtdz_dw * km(i-1)*uv(i-1)+ q*dtdz_up*km(i)*uv(i+1)
!       *(uv(i+1)-uv(i)) !- dt * cmplx(0.0_idx,f) *q1 * uv(i)!&
!            & + dt * uv_adv(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    i=N
    uv(N) = u(N) + cmplx(0.0_idx,1.0_idx) * v(N)
    uv_adv(N) = u_adv(N) + (0.0_idx,1.0_idx) * v_adv(N)
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(i) = cmplx(-1.0_idx * p*dtdz_dw * Km(i-1),0.0_idx)               ! Coefficient of U(i-1)_new
    B(i) = cmplx(1.0_idx + p*dtdz_dw*Km(i-1),p1*f*dt)        ! Coefficient of U(i)_new
    C(i) = cmplx(0.0_idx,0.0_idx)             ! Coefficient of U(i+1)_new
    D(i) = uv(i) + q * dtdz_dw * km(i-1)*(uv(i-1)-uv(i))! &
 !        & +  dt * cmplx(tau_x/rho,tau_y/rho) / (zz(N)-zz(N-1)) &
 !        & - dt * cmplx(0.0_idx,f) *q1 * uv(i)!&
!         & + dt * uv_adv(i)  ! Coefficient of RHS
    ! At the bottom
    i=1
    uv(1) = u(1) + cmplx(0.0_idx,1.0_idx) * v(1)
    uv_adv(1) = u_adv(1) + (0.0_idx,1.0_idx) * v_adv(1)
#ifdef zero_bottom
    A(1) = cmplx(0.0_idx,0.0_idx)  ! Coefficient of U(0),V(0)<-does not exist
    B(1) = cmplx(1.0_idx,0.0_idx) ! Coefficient of U(1),V(1)
    C(1) = cmplx(-1.0_idx,0.0_idx)  ! Coefficient of U(2),V(2)
    D(1) = cmplx(0.0_idx,0.0_idx)  ! Ensure that U(1)=U(2); V(1)=V(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(i) = cmplx(0.0_idx,0.0_idx)               ! Coefficient of U(i-1)_new
    B(i) = cmplx(1.0_idx + p*dtdz_up*Km(i),p1*f*dt)        ! Coefficient of U(i)_new
!    B(i) = cmplx(1.0_idx + p*dtdz_up*Km(i),0.0_idx)        ! Coefficient of U(i)_new
    C(i) = cmplx(-1.0_idx * p* dtdz_up * Km(i),0.0_idx)             ! Coefficient of U(i+1)_new
    D(i) = uv(i) + q * dtdz_up * km(i)*(uv(i+1)-uv(i)) !&
    !     & - dt * cmplx(0.0_idx,f) *q1 * uv(i) !&
    !       & + dt * uv_adv(i)               ! Coefficient of RHS
#endif
    call solve_tri_implicit_complex_nolap(N,A,B,C,D,uv_next)
    u_next=real(uv_next) ; v_next=imag(uv_next)
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))       
       u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
       v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
       u_cor(i) = u_cor(i) + f * v_next(i)
       v_cor(i) = v_cor(i) - f * u_next(i)
       u_vdiff(i) = u_vdiff(i) + Km(i) * dtdz_up  * (p*u_next(i+1)+q*u(i+1)-p*u_next(i)-q*u(i)) - &
            & Km(i-1) * dtdz_dw  * (p*u_next(i)+q*u(i)-p*u_next(i-1)-q*u(i-1))
       v_vdiff(i) = v_vdiff(i) + Km(i) * dtdz_up  * (p*v_next(i+1)+q*v(i+1)-p*v_next(i)-q*v(i)) - &
            & Km(i-1) * dtdz_dw  * (p*v_next(i)+q*v(i)-p*v_next(i-1)-q*v(i-1))
       u_relax(i) = u_relax(i) + u_adv(i)
       v_relax(i) = v_relax(i) + v_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
    v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
    u_cor(i) = u_cor(i) + f * v_next(i)
    v_cor(i) = v_cor(i) - f * u_next(i)
    u_vdiff(i) = u_vdiff(i) + (tau_x/rho) / (zz(N)-zz(N-1)) - &
         & Km(i-1) * dtdz_dw  * (p*u_next(i)+q*u(i)-p*u_next(i-1)-q*u(i-1))
    v_vdiff(i) = v_vdiff(i) + (tau_y/rho) / (zz(N)-zz(N-1)) - &
         & Km(i-1) * dtdz_dw  * (p*v_next(i)+q*v(i)-p*v_next(i-1)-q*v(i-1))
    u_relax(i) = u_relax(i) + u_adv(i)
    v_relax(i) = v_relax(i) + v_adv(i)
    i=1
#ifdef zero_bottom
    u_rate(i) = u_rate(i+1); u_cor(i) = u_cor(i+1)
    u_vdiff(i) = u_vdiff(i+1); u_relax(i) = u_relax(i+1)
    v_rate(i) = v_rate(i+1);  v_cor(i) = v_cor(i+1)
    v_vdiff(i) = v_vdiff(i+1); v_relax(i) = v_relax(i+1)
#else
    dtdz_up = 1.0_idx / ((zz(2)-zz(1))*(depth(i+1)-depth(i)))       
    u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
    v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
    u_cor(i) = u_cor(i) + f * v_next(i)
    v_cor(i) = v_cor(i) - f * u_next(i)
    u_vdiff(i) = u_vdiff(i) + km(i) * dtdz_up  * (p*u_next(i+1)+q*u(i+1)-p*u_next(i)-q*u(i))
    v_vdiff(i) = v_vdiff(i) + km(i) * dtdz_up  * (p*v_next(i+1)+q*v(i+1)-p*v_next(i)-q*v(i))
    u_relax(i) = u_relax(i) + u_adv(i); v_relax(i) = v_relax(i) + v_adv(i)
#endif
  end subroutine cal_next_uv_cn_bud

  subroutine cal_next_uv_dev_bud(dt,N,depth,zz,u,v,km,f,tau_x,tau_y,u_adv,v_adv,u_next,v_next,&
       & u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),u(N),v(N),km(N-1)
    real(idx),intent(in) :: f,tau_x,tau_y
    real(idx),intent(in) :: u_adv(N),v_adv(N)
    real(idx),intent(inout) :: u_next(N),v_next(N)
    real(idx),intent(inout) :: u_rate(n),u_cor(n),u_vdiff(n),u_relax(n)
    real(idx),intent(inout) :: v_rate(n),v_cor(n),v_vdiff(n),v_relax(n)
    real(idx) :: A(N),B(N),C(N),D(N)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i
    real(idx),parameter :: p=1.0_idx,q=1.0_idx-p
    real(idx) :: u_tmp,v_tmp
    u_next = 0.0_idx; v_next=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = -1.0_idx * p*dtdz_dw * Km(i-1)                ! Coefficient of U(i-1)_new
       B(i) = 1.0_idx + p*(dtdz_dw*Km(i-1)+dtdz_up*Km(i))   ! Coefficient of U(i)_new
       C(i) = -1.0_idx * p* dtdz_up * Km(i)                 ! Coefficient of U(i+1)_new
       D(i) = u(i) + q * dtdz_dw * km(i-1)*(u(i-1)-u(i)) &
            & + q * dtdz_up * km(i)*(u(i+1)-u(i)) &
            ! & + dt * 0.5_idx * f * v(i) &
            & + dt * u_adv(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    i=N
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(i) = -1.0_idx * p*dtdz_dw * Km(i-1)                ! Coefficient of U(i-1)_new
    B(i) = 1.0_idx + p*dtdz_dw*Km(i-1)   ! Coefficient of U(i)_new
    C(i) = 0.0_idx
    D(i) = u(i) + q * dtdz_dw * km(i-1)*(u(i-1)-u(i)) &
         ! & + dt * 0.5_idx * f * v(i) &
         & + dt * u_adv(i) &              ! Coefficient of RHS
         & + dt * (tau_x/rho) / (zz(N)-zz(N-1)) 
    ! At the bottom
    i=1
#ifdef zero_bottom
    A(1) = 0.0_idx   ! Coefficient of U(0),V(0)<-does not exist
    B(1) = -1.0_idx  ! Coefficient of U(1),V(1)
    C(1) = 1.0_idx   ! Coefficient of U(2),V(2)
    D(1) = 0.0_idx   ! Ensure that U(1)=U(2); V(1)=V(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(i) = 0.0_idx        ! Coefficient of U(i-1)_new
    B(i) = 1.0_idx + p*dtdz_up*Km(i)   ! Coefficient of U(i)_new
    C(i) = -1.0_idx * p* dtdz_up * Km(i)                 ! Coefficient of U(i+1)_new
    D(i) = u(i) + q * dtdz_up * km(i)*(u(i+1)-u(i)) &
!         & + dt * 0.5_idx * f * v(i)&
         & + dt * u_adv(i)               ! Coefficient of RHS
#endif
    call solve_tri_implicit_real_nolap(N,A,B,C,D,u_next)
    ! V equation
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       A(i) = -1.0_idx * p*dtdz_dw * Km(i-1)                ! Coefficient of U(i-1)_new
       B(i) = 1.0_idx + p*(dtdz_dw*Km(i-1)+dtdz_up*Km(i))   ! Coefficient of U(i)_new
       C(i) = -1.0_idx * p* dtdz_up * Km(i)                 ! Coefficient of U(i+1)_new
       D(i) = u(i) + q * dtdz_dw * km(i-1)*(v(i-1)-v(i)) &
            & + q * dtdz_up * km(i)*(v(i+1)-v(i)) !&
            !            & - 0.5_idx * dt * f * u(i)&
            !& + dt * v_adv(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    i=N
    dtdz_dw = dt / ((zz(N)-zz(N-1))*(depth(N)-depth(N-1)))
    A(i) = -1.0_idx * p*dtdz_dw * Km(i-1)                ! Coefficient of U(i-1)_new
    B(i) = 1.0_idx + p*dtdz_dw*Km(i-1)   ! Coefficient of U(i)_new
    C(i) = 0.0_idx
    D(i) = v(i) + q * dtdz_dw * km(i-1)*(v(i-1)-v(i)) !&
         !         & - dt * 0.5_idx * f * u(i) &
         !        & + dt * v_adv(i) &              ! Coefficient of RHS
         !        & +  dt * (tau_y/rho) / (zz(N)-zz(N-1)) 
    ! At the bottom
    i=1
#ifdef zero_bottom
    A(1) = 0.0_idx   ! Coefficient of U(0),V(0)<-does not exist
    B(1) = -1.0_idx  ! Coefficient of U(1),V(1)
    C(1) = 1.0_idx   ! Coefficient of U(2),V(2)
    D(1) = 0.0_idx   ! Ensure that U(1)=U(2); V(1)=V(2)
#else
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(i) = 0.0_idx        ! Coefficient of U(i-1)_new
    B(i) = 1.0_idx + p * dtdz_up*Km(i)   ! Coefficient of U(i)_new
    C(i) = -1.0_idx * p * dtdz_up * Km(i)                 ! Coefficient of U(i+1)_new
    D(i) = v(i) + q * dtdz_up * km(i)*(v(i+1)-v(i)) &
         !& - dt * 0.5_idx * f * u(i)&
         & + dt * v_adv(i)               ! Coefficient of RHS
#endif
    call solve_tri_implicit_real_nolap(N,A,B,C,D,v_next)
!    do i=2,N
!       u_tmp=u_next(i); v_tmp=v_next(i)
!       u_next(i)=u_next(i) + dt * 0.5_idx*f*v_tmp
!       v_next(i)=v_next(i) - dt * 0.5_idx*f*u_tmp
!    end do

    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))       
       u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
       v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
       u_cor(i) = u_cor(i) + f * v_next(i)
       v_cor(i) = v_cor(i) - f * u_next(i)
       u_vdiff(i) = u_vdiff(i) + Km(i) * dtdz_up  * (p*u_next(i+1)+q*u(i+1)-p*u_next(i)-q*u(i)) - &
            & Km(i-1) * dtdz_dw  * (p*u_next(i)+q*u(i)-p*u_next(i-1)-q*u(i-1))
       v_vdiff(i) = v_vdiff(i) + Km(i) * dtdz_up  * (p*v_next(i+1)+q*v(i+1)-p*v_next(i)-q*v(i)) - &
            & Km(i-1) * dtdz_dw  * (p*v_next(i)+q*v(i)-p*v_next(i-1)-q*v(i-1))
       u_relax(i) = u_relax(i) + u_adv(i)
       v_relax(i) = v_relax(i) + v_adv(i)
    end do
    i=N
    dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
    u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
    v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
    u_cor(i) = u_cor(i) + f * v_next(i)
    v_cor(i) = v_cor(i) - f * u_next(i)
    u_vdiff(i) = u_vdiff(i) + (tau_x/rho) / (zz(N)-zz(N-1)) - &
         & Km(i-1) * dtdz_dw  * (p*u_next(i)+q*u(i)-p*u_next(i-1)-q*u(i-1))
    v_vdiff(i) = v_vdiff(i) + (tau_y/rho) / (zz(N)-zz(N-1)) - &
         & Km(i-1) * dtdz_dw  * (p*v_next(i)+q*v(i)-p*v_next(i-1)-q*v(i-1))
    u_relax(i) = u_relax(i) + u_adv(i)
    v_relax(i) = v_relax(i) + v_adv(i)
    i=1
#ifdef zero_bottom
    u_rate(i) = u_rate(i+1); u_cor(i) = u_cor(i+1)
    u_vdiff(i) = u_vdiff(i+1); u_relax(i) = u_relax(i+1)
    v_rate(i) = v_rate(i+1);  v_cor(i) = v_cor(i+1)
    v_vdiff(i) = v_vdiff(i+1); v_relax(i) = v_relax(i+1)
#else
    dtdz_up = 1.0_idx / ((zz(2)-zz(1))*(depth(i+1)-depth(i)))       
    u_rate(i) = u_rate(i) + (u_next(i) - u(i)) / dt
    v_rate(i) = v_rate(i) + (v_next(i) - v(i)) / dt
    u_cor(i) = u_cor(i) + f * v_next(i)
    v_cor(i) = v_cor(i) - f * u_next(i)
    u_vdiff(i) = u_vdiff(i) + km(i) * dtdz_up  * (p*u_next(i+1)+q*u(i+1)-p*u_next(i)-q*u(i))
    v_vdiff(i) = v_vdiff(i) + km(i) * dtdz_up  * (p*v_next(i+1)+q*v(i+1)-p*v_next(i)-q*v(i))
    u_relax(i) = u_relax(i) + u_adv(i); v_relax(i) = v_relax(i) + v_adv(i)
#endif
  end subroutine cal_next_uv_dev_bud

  !---------------------------------
  ! TKE
  !---------------------------------
  ! subroutine for calculating new tke
  subroutine cal_next_qq(dt,N,depth,zz,qq,bvf,shear2,L,kq,Km,Kt,tau_x,tau_y,B_1,qq_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),qq(N),kq(N-1),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N),tau_x,tau_y
    real(idx),intent(in) :: B_1
    real(idx),intent(inout) :: qq_next(N)
    real(idx) :: RHS_qq(N-1),SP(N-1),BP(N-1)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up,tau_ts
    integer :: i
    RHS_qq = 0.0_idx ; SP = 0.0_idx ; BP=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx

    ! Make matrix
    do i = 2,N-1
       dtdz_dw = dt / ((depth(i+1)-depth(i))*(zz(i)-zz(i-1)))
       dtdz_up = dt / ((depth(i+1)-depth(i))*(zz(i+1)-zz(i)))
       A(i) = - 1.0_idx * dtdz_dw * kq(i-1)                     ! Coefficient of qq(i-1)_new
       B(i) = 1.0_idx + (dtdz_dw * kq(i-1) + dtdz_up * kq(i))   ! Coefficient of qq(i)_new
       C(i) = - 1.0_idx * (dtdz_up * kq(i))                ! Coefficient of qq(i+1)_new
       !------------------------------------------------
       ! Compute Shear production and Buoyancy production
       !------------------------------------------------
       ! Buoyancy frequency
       SP(i) = km(i) * Shear2(i)
       ! ROMS correction
       !       if (bvf(i) .lt. 0.0_idx  .and. bvf(i) .ge. -5.0e-5_idx) then
       !          BP(i) = 0.0_idx
       !       else          
       BP(i) = -1.0_idx * kt(i) * bvf(i)
       !      end if
       RHS_qq(i) = dt * (2.0_idx * (SP(i)+BP(i)) - (2.0_idx * qq(i)**1.5_idx) / (B_1 * l(i))) ! Forcing term
       D(i) = qq(i) +  RHS_qq(i)                  ! Coefficient of RHS
    end do
    ! Boundary Condition---------------------------------
    ! At the surface-----------
    tau_ts=sqrt(tau_x**2+tau_y**2) / rho ! = u_star^2 (tau / rho)
    A(N) = 0.0_idx ! Coefficient of qq_new(N-1) 
    B(N) = 1.0_idx ! Coefficient of qq_new(i)
    C(N) = 0.0_idx ! Coefficient of qq_new(N+1) <-does not exist
    D(N) = ((B_1)**(2.0_idx/3.0_idx))*(tau_ts)! surface energy
    ! At the bottom------------
    A(1) = 0.0_idx   ! Coefficient of qq(0)_new
    B(1) = 1.0_idx   ! Coefficient of qq(1)_new
    C(1) = 0.0_idx   ! Coefficient of qq(2)_new
    D(1) = 0.0_idx
    !----------------------------------------------------
    ! Solve---------------------------------------------
    call solve_tri_implicit_real_nolap(N,A,B,C,D,qq_next)
    do i=1,N
       if(qq_next(i) .le. qq_min) then
          qq_next(i) = qq_min
       end if
    end do
  end subroutine cal_next_qq
  subroutine cal_next_qq_bud(dt,N,depth,zz,qq,bvf,shear2,L,kq,Km,Kt,&
       & tau_x,tau_y,B_1,qq_next,qq_tend,qq_vdiff,qq_sp,qq_bp,qq_disp)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),qq(N),kq(N-1),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N),tau_x,tau_y
    real(idx),intent(in) :: B_1
    real(idx),intent(inout) :: qq_next(N)
    real(idx),intent(inout) :: qq_tend(N),qq_vdiff(N),qq_sp(N),qq_bp(N),qq_disp(N)    
    real(idx) :: RHS_qq(N-1),SP(N-1),BP(N-1)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up,tau_ts
    integer :: i
    RHS_qq = 0.0_idx ; SP = 0.0_idx ; BP=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       dtdz_dw = dt / ((depth(i+1)-depth(i))*(zz(i)-zz(i-1)))
       dtdz_up = dt / ((depth(i+1)-depth(i))*(zz(i+1)-zz(i)))
       A(i) = - 1.0_idx * dtdz_dw * kq(i-1)                     ! Coefficient of qq(i-1)_new
       B(i) = 1.0_idx + (dtdz_dw * kq(i-1) + dtdz_up * kq(i))   ! Coefficient of qq(i)_new
       C(i) = - 1.0_idx * (dtdz_up * kq(i))                ! Coefficient of qq(i+1)_new
       !------------------------------------------------
       ! Compute Shear production and Buoyancy production
       !------------------------------------------------
       ! Buoyancy frequency
       SP(i) = km(i) * Shear2(i)
       BP(i) = -1.0_idx * kt(i) * bvf(i)
       RHS_qq(i) = dt * (2.0_idx * (SP(i)+BP(i)) - (2.0_idx * qq(i)**1.5_idx) / (B_1 * l(i))) ! Forcing term
       D(i) = qq(i) +  RHS_qq(i)                  ! Coefficient of RHS
    end do
    ! Boundary Condition---------------------------------
    ! At the surface-----------
    tau_ts=sqrt(tau_x**2+tau_y**2) / rho ! = u_star^2 (tau / rho)
    A(N) = 0.0_idx ! Coefficient of qq_new(N-1) 
    B(N) = 1.0_idx ! Coefficient of qq_new(i)
    C(N) = 0.0_idx ! Coefficient of qq_new(N+1) <-does not exist
    D(N) = ((B_1)**(2.0_idx/3.0_idx))*(tau_ts)! surface energy
    ! At the bottom------------
    A(1) = 0.0_idx   ! Coefficient of qq(0)_new
    B(1) = 1.0_idx   ! Coefficient of qq(1)_new
    C(1) = 0.0_idx   ! Coefficient of qq(2)_new
    D(1) = 0.0_idx
    !----------------------------------------------------
    ! Solve---------------------------------------------
    call solve_tri_implicit_real_nolap(N,A,B,C,D,qq_next)
    do i=1,N
       if(qq_next(i) .le. qq_min) then
          qq_next(i) = qq_min
       end if
    end do
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       qq_tend(i) = qq_tend(i) + (qq_next(i) - qq(i)) / dt
       qq_vdiff(i) = qq_vdiff(i) + Kq(i) * dtdz_up  * (QQ_next(i+1)-QQ_next(i)) - &
            & Kq(i-1) * dtdz_dw  * (QQ_next(i)-QQ_next(i-1))
       qq_sp(i)=qq_sp(i)+2.0_idx * SP(i)
       qq_bp(i)=qq_bp(i)+2.0_idx * BP(i)
       qq_disp(i) = qq_disp(i) - 2.0_idx * (qq(i)**1.5_idx) / (B_1 * l(i))
    end do
  end subroutine cal_next_qq_bud

  subroutine cal_next_qq_cn_bud(dt,N,depth,zz,qq,bvf,shear2,L,kq,Km,Kt,&
       & tau_x,tau_y,B_1,qq_next,qq_tend,qq_vdiff,qq_sp,qq_bp,qq_disp)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),qq(N),kq(N-1),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N),tau_x,tau_y
    real(idx),intent(in) :: B_1
    real(idx),intent(inout) :: qq_next(N)
    real(idx),intent(inout) :: qq_tend(N),qq_vdiff(N),qq_sp(N),qq_bp(N),qq_disp(N)    
    real(idx) :: RHS_qq(N-1),SP(N-1),BP(N-1)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up,tau_ts
    integer :: i
    real(idx),parameter :: p=0.5_idx,q=1.0_idx-p
    RHS_qq = 0.0_idx ; SP = 0.0_idx ; BP=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       dtdz_dw = dt / ((depth(i+1)-depth(i))*(zz(i)-zz(i-1)))
       dtdz_up = dt / ((depth(i+1)-depth(i))*(zz(i+1)-zz(i)))
       A(i) = - 1.0_idx * p * dtdz_dw * kq(i-1)                     ! Coefficient of qq(i-1)_new
       B(i) = 1.0_idx + p *(dtdz_dw * kq(i-1) + dtdz_up * kq(i))   ! Coefficient of qq(i)_new
       C(i) = - 1.0_idx * p * dtdz_up * kq(i)                ! Coefficient of qq(i+1)_new
       !------------------------------------------------
       ! Compute Shear production and Buoyancy production
       !------------------------------------------------
       ! Buoyancy frequency
       SP(i) = km(i) * Shear2(i)
       BP(i) = -1.0_idx * kt(i) * bvf(i)
       RHS_qq(i) = dt * (2.0_idx * (SP(i)+BP(i)) - (2.0_idx * qq(i)**1.5_idx) / (B_1 * l(i))) ! Forcing term
       D(i) = qq(i) +  RHS_qq(i) + q*dtdz_dw * kq(i-1)*(qq(i-1)-qq(i))+ q*dtdz_up * kq(i)*(qq(i+1)-qq(i))                  ! Coefficient of RHS
    end do
    ! Boundary Condition---------------------------------
    ! At the surface-----------
    tau_ts=sqrt(tau_x**2+tau_y**2) / rho ! = u_star^2 (tau / rho)
    A(N) = 0.0_idx ! Coefficient of qq_new(N-1) 
    B(N) = 1.0_idx ! Coefficient of qq_new(i)
    C(N) = 0.0_idx ! Coefficient of qq_new(N+1) <-does not exist
    D(N) = ((B_1)**(2.0_idx/3.0_idx))*(tau_ts)! surface energy
    ! At the bottom------------
    A(1) = 0.0_idx   ! Coefficient of qq(0)_new
    B(1) = 1.0_idx   ! Coefficient of qq(1)_new
    C(1) = 0.0_idx   ! Coefficient of qq(2)_new
    D(1) = 0.0_idx
    !----------------------------------------------------
    ! Solve---------------------------------------------
    call solve_tri_implicit_real_nolap(N,A,B,C,D,qq_next)
    do i=1,N
       if(qq_next(i) .le. qq_min) then
          qq_next(i) = qq_min
       end if
    end do
    ! Budget===========================
    do i = 2,N-1
       dtdz_dw = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i)-depth(i-1)))
       dtdz_up = 1.0_idx / ((zz(i)-zz(i-1))*(depth(i+1)-depth(i)))
       qq_tend(i) = qq_tend(i) + (qq_next(i) - qq(i)) / dt
       qq_vdiff(i) = qq_vdiff(i) + Kq(i) * dtdz_up  * (p*qq_next(i+1)+q*qq(i+1)-p*qq_next(i)-q*qq(i)) - &
            & Kq(i-1) * dtdz_dw  * (p*qq_next(i)+q*qq(i)-p*qq_next(i-1)-q*qq(i-1))
       qq_sp(i)=qq_sp(i)+2.0_idx * SP(i)
       qq_bp(i)=qq_bp(i)+2.0_idx * BP(i)
       qq_disp(i) = qq_disp(i) - 2.0_idx * (qq(i)**1.5_idx) / (B_1 * l(i))
    end do
  end subroutine cal_next_qq_cn_bud

  ! Calculate next length scale
  subroutine cal_next_l2(dt,N,depth,zz,qq,bvf,shear2,L,kq,Km,Kt,qq_next,B_1,E_1,E_2,l_next)
    implicit none
    integer,parameter :: idx=8
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),qq(N),kq(N-1),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N)
    real(idx),intent(in) :: qq_next(N)
    real(idx),intent(in) :: B_1,E_1,E_2
    real(idx),intent(inout) :: l_next(N)
    real(idx) :: RHS_qql(N-1),SP(N-1),BP(N-1),qql_next(N)
    real(idx) :: A(N),B(N),C(N),D(N),dtdz_dw,dtdz_up
    real(idx),parameter :: kappa=0.41_idx
    integer :: i
    RHS_qql = 0.0_idx ; SP = 0.0_idx ; BP=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx

    ! Make matrix
    do i = 2,N-1
       dtdz_dw = dt / ((depth(i+1)-depth(i))*(zz(i)-zz(i-1)))
       dtdz_up = dt / ((depth(i+1)-depth(i))*(zz(i+1)-zz(i)))
       A(i) = - 1.0_idx * dtdz_dw * kq(i-1)                     ! Coefficient of qq(i-1)_new
       B(i) = 1.0_idx + (dtdz_dw * kq(i-1) + dtdz_up * kq(i))   ! Coefficient of qq(i)_new
       C(i) = - 1.0_idx * (dtdz_up * kq(i))                ! Coefficient of qq(i+1)_new
       !------------------------------------------------
       ! Compute Shear production and Buoyancy production
       !------------------------------------------------
       ! Buoyancy frequency
       SP(i) = km(i) * Shear2(i)
       BP(i) = -1.0_idx * kt(i) * bvf(i)
       RHS_qql(i) = dt * (L(i) * E_1 * (SP(i)+BP(i)) - (qq(i)**1.5_idx / B_1) * &
            & (1.0_idx + E_2 * (l(i)  / (kappa * abs(zz(i))))**2))
       D(i) = qq(i) +  RHS_qql(i)                  ! Coefficient of RHS
    end do
    ! Boundary Condition---------------------------------
    ! At the surface-----------
    A(N) = 0.0_idx ! Coefficient of qq_new(N-1) 
    B(N) = 1.0_idx ! Coefficient of qq_new(i)
    C(N) = 0.0_idx ! Coefficient of qq_new(N+1) <-does not exist
    D(N) = 0.0_idx
    ! At the bottom------------
    A(1) = 0.0_idx   ! Coefficient of qq(0)_new
    B(1) = 1.0_idx   ! Coefficient of qq(1)_new
    C(1) = 0.0_idx   ! Coefficient of qq(2)_new
    D(1) = 0.0_idx
    !----------------------------------------------------
    ! Solve---------------------------------------------
    call solve_tri_implicit_real_nolap(N,A,B,C,D,qql_next)
    do i=1,N
       if(qql_next(i) .le. qql_min) then
          qql_next(i) = qql_min
       end if
       l_next(i)= qql_next(i) / qq_next(i)
    end do
  end subroutine cal_next_l2
  subroutine cal_next_l(dt,N,depth,zz,qq,bvf,shear2,L,kq,Km,Kt,qq_next,B_1,E_1,E_2,l_next)
    implicit none
    integer,parameter :: idx=8
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N),qq(N),kq(N-1),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N)
    real(idx),intent(in) :: qq_next(N)
    real(idx),intent(in) :: B_1,E_1,E_2
    real(idx),intent(inout) :: l_next(N)
    real(idx),parameter :: kappa=0.41_idx
    integer :: iz
    real(idx) :: qz_int,q_int,lt_inv
    do iz = n-1,1,-1
       ! Integration
       q_int=q_int + (zz(iz+1)-zz(iz))*(sqrt(qq(iz+1))+sqrt(qq(iz)))*0.5_idx
       qz_int=qz_int + (zz(iz+1)-zz(iz))*(sqrt(qq(iz+1))*abs(zz(iz+1))+sqrt(qq(iz))*abs(zz(iz)))*0.5
    end do
    lt_inv = 0.1_idx * qz_int / q_int
    do iz = 1,n
       l_next(iz)=  lt_inv * kappa * abs(zz(iz)) / (lt_inv+kappa * abs(zz(iz)))
    end do
  end subroutine cal_next_l
end module solve_diag
