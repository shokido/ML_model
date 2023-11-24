module solve_diag
  use ml_param
  implicit none
  private
  !-------------------------------------------------------
  public :: solve_tri_implicit_real_nolap, solve_tri_implicit_complex_nolap
  public :: cal_next_theta
  public :: cal_next_salt
  public :: cal_next_uv
  public :: cal_next_qq
contains
  !============================================================================
  subroutine solve_tri_implicit_real_nolap(N,A_in,B_in,C_in,D_in,U)
    ! Solve A * U_new(i-1) + B * U_new(i) + C * U_new(i+1) = D
    ! solve DT/dt=(k(z)(d^2 T)/(dz^2))
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: A_in(1:N),B_in(1:N),C_in(1:N),D_in(1:N)
    real(idx),intent(inout) :: U(1:N)
    real(idx),allocatable :: A(:),B(:),C(:),D(:)
    integer :: i
    real(idx) :: m
    allocate(A(1:N));allocate(B(1:N));allocate(C(1:N));allocate(D(1:N))
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
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
  end subroutine solve_tri_implicit_real_nolap
  !============================================================================
  subroutine solve_tri_implicit_complex_nolap(N,A_in,B_in,C_in,D_in,U)
    ! Solve A * U_new(i-1) + B * U_new(i) + C * U_new(i+1) = D
    ! solve DT/dt=(k(z)(d^2 T)/(dz^2))
    implicit none
    integer,intent(in) :: N
    complex(kind(0d0)),intent(in) :: A_in(1:N),B_in(1:N),C_in(1:N),D_in(1:N)
    complex(kind(0d0)),intent(inout) :: U(N)
    complex(kind(0d0)),allocatable :: A(:),B(:),C(:),D(:)
    integer :: i
    complex(kind(0d0)) :: m
    allocate(A(1:N));allocate(B(1:N));allocate(C(1:N));allocate(D(1:N))
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
  end subroutine solve_tri_implicit_complex_nolap

  !============================================================================
  ! subroutine for calculating new temperature
  subroutine cal_next_theta(dt,N,z_rho,z_q,theta,kt,heat,heat_no_solar,heat_correct,theta_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,z_rho(N),z_q(N+1),theta(N),kt(N-1),heat(N+1),heat_correct(N)
    real(idx),intent(inout) :: theta_next(N)
    real(idx),intent(in) :: heat_no_solar
    integer :: i
    integer :: NB
    real(idx) :: dtdz_dw,dtdz_up
    real(idx),allocatable :: A(:),B(:),C(:),D(:),tmp(:)
    NB=N
    allocate(A(1:NB));allocate(B(1:NB));allocate(C(1:NB));allocate(D(1:NB))
    allocate(tmp(1:NB))
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    tmp=0.0_idx
    ! Make matrix------------------------
    do i = 2,NB-1
       dtdz_up = dt / ((z_q(i)-z_q(i+1))*(z_rho(i-1)-z_rho(i)))
       dtdz_dw = dt / ((z_q(i)-z_q(i+1))*(z_rho(i)-z_rho(i+1)))
       A(i) = - 1.0_idx * dtdz_up * Kt(i-1)      ! Coeff. of T(i-1)_new
       B(i) = 1.0_idx + (dtdz_dw * Kt(i)+dtdz_up*Kt(i-1))    ! Coefficient of T(i)_new
       C(i) = - 1.0_idx * dtdz_dw * Kt(i)                 ! Coefficient of T(i+1)_new
       D(i) = theta(i) + dt * (heat(i-1)-heat(i)) / (rho*cpw*(z_q(i)-z_q(i+1)))  &
            & +dt * heat_correct(i)  ! Coefficient of RHS
    end do
    ! Boundary Condition----------------
    ! At the surface
    ! KT(N)*dT/dz = Q
    dtdz_dw = dt / ((z_q(1)-z_q(2))*(z_rho(1)-z_rho(2)))
    C(1) = -1.0_idx * dtdz_dw * Kt(1)
    B(1) = 1.0_idx -C(1)
    A(1) = 0.0_idx
    D(1) =  theta(1) + dt * heat_no_solar / (rho *cpw*(z_q(1)-z_q(2))) + &
         & dt * (heat(1)-heat(2)) / (rho*cpw*(z_q(1)-z_q(2)))+ dt * heat_correct(1)
    ! At the bottom
    dtdz_up = dt / ((z_q(NB)-z_q(NB+1))*(z_rho(NB-1)-z_rho(NB)))
    A(NB) = - 1.0_idx * dtdz_up * Kt(1)
    B(NB) = 1.0_idx-A(NB)
    C(NB) = 0.0_idx
    D(NB) = theta(NB)+dt*heat_correct(NB)   ! Coefficient of RHS
    call solve_tri_implicit_real_nolap(NB,A,B,C,D,tmp)
    theta_next(1:NB)=tmp(1:NB)
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
  end subroutine cal_next_theta
  !---------------------------------
  ! Salinity
  !---------------------------------
  subroutine cal_next_salt(dt,N,z_rho,z_q,salt,ks,ssflux,salt_correct,salt_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,z_rho(N),z_q(N+1),salt(N),ks(N-1),salt_correct(N),ssflux
    real(idx),intent(inout) :: salt_next(N)
    integer :: i
    integer :: NB
    real(idx) :: dtdz_dw,dtdz_up
    real(idx),allocatable :: A(:),B(:),C(:),D(:),tmp(:)
    NB=N
    allocate(A(1:NB));allocate(B(1:NB));allocate(C(1:NB));allocate(D(1:NB))
    allocate(tmp(1:NB))
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    tmp=0.0_idx
    ! Make matrix------------------------
    do i = 2,NB-1
       dtdz_up = dt / ((z_q(i)-z_q(i+1))*(z_rho(i-1)-z_rho(i)))
       dtdz_dw = dt / ((z_q(i)-z_q(i+1))*(z_rho(i)-z_rho(i+1)))
       A(i) = - 1.0_idx * dtdz_up * Ks(i-1)      ! Coeff. of T(i-1)_new
       B(i) = 1.0_idx + (dtdz_dw * Ks(i)+dtdz_up*Ks(i-1))    ! Coefficient of T(i)_new
       C(i) = - 1.0_idx * dtdz_dw * Ks(i)                 ! Coefficient of T(i+1)_new
       D(i) = salt(i)+dt * salt_correct(i)  ! Coefficient of RHS
    end do
    ! Boundary Condition----------------
    ! At the surface
    ! KS(N)*dT/dz = Q
    dtdz_dw = dt / ((z_q(1)-z_q(2))*(z_rho(1)-z_rho(2)))
    C(1) = -1.0_idx * dtdz_dw * Ks(1)
    B(1) = 1.0_idx - C(1)
    A(1) = 0.0_idx
    D(1) =  salt(1) + dt * ssflux / (z_q(1)-z_q(2)) + dt * salt_correct(1)
    ! At the bottom
    dtdz_up = dt / ((z_q(NB)-z_q(NB+1))*(z_rho(NB-1)-z_rho(NB)))
    A(NB) = - 1.0_idx * dtdz_up * Ks(1)
    B(NB) = 1.0_idx -A(NB)
    C(NB) = 0.0_idx
    D(NB) = salt(NB)+dt*salt_correct(NB)   ! Coefficient of RHS
    call solve_tri_implicit_real_nolap(NB,A,B,C,D,tmp)
    salt_next(1:NB)=tmp(1:NB)
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
  end subroutine cal_next_salt
  !---------------------------------! 
  ! Velocity                        !
  !---------------------------------!
  ! subroutine for calculating new velocity
  subroutine cal_next_uv(dt,N,z_rho,z_q,u,v,km,f,tau_x,tau_y,u_correct,v_correct,u_next,v_next)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: dt,z_rho(n),z_q(n+1),u(n),v(n),km(n-1)
    real(idx),intent(in) :: f,tau_x,tau_y
    real(idx),intent(in) :: u_correct(n),v_correct(n)
    real(idx),intent(inout) :: u_next(n),v_next(n)
    complex(kind(0d0)) :: uv(n),uv_correct(n)
    complex(kind(0d0)),allocatable ::A(:),B(:),C(:),D(:),uv_next(:)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i,NB
    NB=N
    allocate(A(1:NB));allocate(B(1:NB));allocate(C(1:NB));allocate(D(1:NB))
    allocate(uv_next(1:NB))
    uv=0.0_idx ; uv_next = 0.0_idx; uv_correct=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,NB-1
       uv(i) = u(i) + (0.0_idx,1.0_idx) * v(i)
       uv_correct(i) = u_correct(i) + (0.0_idx,1.0_idx) * v_correct(i)
       dtdz_up = dt / ((z_q(i)-z_q(i+1))*(z_rho(i-1)-z_rho(i)))
       dtdz_dw = dt / ((z_q(i)-z_q(i+1))*(z_rho(i)-z_rho(i+1)))
       A(i) = cmplx(-1.0_idx * dtdz_up * Km(i-1),0.0_idx)               ! Coefficient of U(i-1)_new
       B(i) = cmplx(1.0_idx + (dtdz_up*Km(i-1)+dtdz_dw*Km(i)),f*dt)        ! Coefficient of U(i)_new
       C(i) = cmplx(-1.0_idx * dtdz_dw * Km(i),0.0_idx)             ! Coefficient of U(i+1)_new
       D(i) = uv(i)+ dt * uv_correct(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    uv(1) = u(1) + (0.0_idx,1.0_idx) * v(1)
    uv_correct(1) = u_correct(1) + (0.0_idx,1.0_idx) * v_correct(1)
    dtdz_dw = dt / ((z_q(1)-z_q(2))*(z_rho(1)-z_rho(2)))
    A(1) = cmplx(0.0_idx,0_idx)
    C(1) = cmplx(-1.0_idx * dtdz_dw * Km(1),0.0_idx)   ! Coefficient of U,V(N-1)        
    B(1) = cmplx(1.0_idx + dtdz_dw * Km(1),f*dt)       ! Coefficient of U,V(N)
    D(1) = uv(1) +  dt * cmplx(tau_x/rho,tau_y/rho) / (z_q(1)-z_q(2))+ dt * uv_correct(1)
    ! At the bottom
    dtdz_up = dt / ((z_q(NB)-z_q(NB+1))*(z_rho(NB-1)-z_rho(NB)))
    A(NB) = cmplx(-1.0_idx * dtdz_up * Km(NB-1),0.0_idx)
    C(NB) = cmplx(0.0_idx,0.0_idx)
    B(NB) = cmplx(1.0_idx + dtdz_up*Km(NB-1),f*dt)
    D(NB) = uv(NB)+ dt * uv_correct(NB)
    call solve_tri_implicit_complex_nolap(NB,A,B,C,D,uv_next)
    u_next(1:NB)=real(uv_next(1:NB),idx)
    v_next(1:NB)=imag(uv_next(1:NB))
  end subroutine cal_next_uv

  !---------------------------------
  ! TKE
  !---------------------------------
  ! subroutine for calculating new tke
  subroutine cal_next_qq(dt,N,z_rho,z_q,qq,bvf,shear2,L,kq,Km,Kt,tau_x,tau_y,B_1,qq_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,z_rho(N),z_q(N+1),qq(N+1),kq(N),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N+1),tau_x,tau_y
    real(idx),intent(in) :: B_1
    real(idx),intent(inout) :: qq_next(N+1)
    real(idx),allocatable :: A(:),B(:),C(:),D(:)
    real(idx),allocatable :: RHS_qq(:),SP(:),BP(:),DS(:),tmp(:)
    integer :: i,NB
    real(idx) :: dtdz_dw,dtdz_up,tau_ts
    NB=N
    allocate(A(1:NB+1));allocate(B(1:NB+1));allocate(C(1:NB+1));allocate(D(1:NB+1))
    allocate(RHS_qq(1:NB+1));allocate(SP(NB+1));allocate(BP(NB+1));allocate(DS(NB+1))
    allocate(tmp(NB+1))
    RHS_qq = 0.0_idx ; SP = 0.0_idx ; BP=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx

    ! qq(1)       qq(2)........qq(N)...qq(N+1) 
    !       T(1)                     T(N)
    !             kt(1)        kt(N-1)
    !       kq(1)                    kq(N)
    !             bvf(1)       bvf(n-1)
    ! Make matrix
    do i = 2,NB
       dtdz_up = dt / ((z_rho(i-1)-z_rho(i))*(z_q(i-1)-z_q(i)))
       dtdz_dw = dt / ((z_rho(i-1)-z_rho(i))*(z_q(i)-z_q(i+1)))
       A(i) = - 1.0_idx * dtdz_up * kq(i-1)                     ! Coefficient of qq(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * kq(i-1) + dtdz_dw * kq(i))   ! Coefficient of qq(i)_new
       C(i) = - 1.0_idx * (dtdz_dw * kq(i))                ! Coefficient of qq(i+1)_new
       !------------------------------------------------
       ! Compute Shear production and Buoyancy production
       !------------------------------------------------
       ! Buoyancy frequency
       SP(i) = km(i-1) * shear2(i-1)
       BP(i) = -1.0_idx * kt(i-1) * bvf(i-1)
       DS(i) = (qq(i)**1.5_idx) / (B_1 * l(i))
       RHS_qq(i) = dt * 2.0_idx * (SP(i)+BP(i)-DS(i)) ! Forcing term
       D(i) = qq(i) +  RHS_qq(i)                  ! Coefficient of RHS
    end do
    ! Boundary Condition---------------------------------
    ! At the surface-----------
    tau_ts=sqrt(tau_x**2+tau_y**2) / rho ! = u_star^2 (tau / rho)
    A(1) = 0.0_idx ! Coefficient of qq_new(N) 
    B(1) = 1.0_idx ! Coefficient of qq_new(N+1)
    C(1) = 0.0_idx ! Coefficient of qq_new(N+2) <-does not exist
    D(1) = ((B_1)**(2.0_idx/3.0_idx))*(tau_ts)! surface energy
    ! At the bottom------------
    A(NB+1) = 0.0_idx   ! Coefficient of qq(0)_new
    B(NB+1) = 1.0_idx   ! Coefficient of qq(1)_new
    C(NB+1) = 0.0_idx   ! Coefficient of qq(2)_new
    D(NB+1) = 0.0_idx
    !----------------------------------------------------
    ! Solve---------------------------------------------
    call solve_tri_implicit_real_nolap(NB+1,A,B,C,D,tmp)
    do i=1,NB+1
       if(tmp(i) .le. qq_min) then
          tmp(i) = qq_min
       end if
    end do
    qq_next(1:NB)=tmp(1:NB)
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
    deallocate(RHS_qq);deallocate(SP);deallocate(BP);deallocate(DS);deallocate(tmp)
  end subroutine cal_next_qq
end module solve_diag
