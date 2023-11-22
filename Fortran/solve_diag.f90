module solve_diag
  use param
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
  subroutine cal_next_theta(dt,N,depth,zz,theta,kt,heat,heat_no_solar,heat_correct,theta_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N+1),theta(N),kt(N-1),heat(N+1),heat_correct(N)
    real(idx),intent(inout) :: theta_next(N)
    real(idx),intent(in) :: heat_no_solar
    integer :: i
    real(idx) :: dtdz_dw,dtdz_up
    real(idx),allocatable :: A(:),B(:),C(:),D(:)
    allocate(A(1:N));allocate(B(1:N));allocate(C(1:N));allocate(D(1:N))
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx

    ! zz(i)---z(i)--zz(i+1)---z(i+1)--
    ! kt(i-1)        kt(i)          
    ! Make matrix------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i+1)-zz(i))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i+1)-zz(i))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * dtdz_dw * Kt(i-1)      ! Coeff. of T(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * Kt(i)+dtdz_dw*Kt(i-1))    ! Coefficient of T(i)_new
       C(i) = - 1.0_idx * dtdz_up * Kt(i)                 ! Coefficient of T(i+1)_new
       D(i) = Theta(i) + dt * (heat(i)-heat(i-1)) / (rho*cpw*(zz(i+1)-zz(i)))  &
            & +dt * heat_correct(i)  ! Coefficient of RHS
    end do
    ! Boundary Condition----------------
    ! At the surface
    dtdz_dw = dt / ((zz(N+1)-zz(N))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * dtdz_dw * Kt(N-1)       ! Coefficient of T(N-1)
    B(N) = 1.0_idx + dtdz_dw * Kt(N-1)        ! Coefficient of T(N)
    C(N) = 0.0_idx                         ! Coefficient of T(N+1)_new(does not exist)
    ! KT(N)*dT/dz = Q
    ! dT/dt = ((Q/rcph)-vdiff(N-1))+() /tau_rst
    D(N) =  Theta(N) + dt * heat_no_solar / (rho *cpw*(zz(N+1)-zz(N))) + &
         & dt * (heat(N)-heat(N-1)) / (rho*cpw*(zz(N+1)-zz(N)))+ dt * heat_correct(N)
    ! At the bottom
    dtdz_up = dt / ((zz(2)-zz(1))*(depth(2)-depth(1)))
    A(1) = 0.0_idx      ! Coeff. of T(0)_new
    B(1) = 1.0_idx + dtdz_up * Kt(1)    ! Coefficient of T(i)_new
    C(1) = - 1.0_idx * dtdz_up * Kt(1)                 ! Coefficient of T(i+1)_new
    D(1) = Theta(1)   ! Coefficient of RHS
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Theta_next)
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
  end subroutine cal_next_theta
  !---------------------------------
  ! Salinity
  !---------------------------------
  ! subroutine for calculating new salinity
  subroutine cal_next_salt(dt,N,depth,zz,salt,ks,ssflux,salt_correct,salt_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N+1),salt(N),ks(N-1),salt_correct(N),ssflux
    real(idx),intent(inout) :: salt_next(N)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i
    real(idx),allocatable :: A(:),B(:),C(:),D(:)
    allocate(A(1:N));allocate(B(1:N));allocate(C(1:N));allocate(D(1:N))
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix-------------------------------------------
    do i = 2,N-1
       dtdz_dw = dt / ((zz(i+1)-zz(i))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i+1)-zz(i))*(depth(i+1)-depth(i)))
       A(i) = - 1.0_idx * dtdz_dw * ks(i-1)               ! Coefficient of S(i-1)_new
       B(i) = 1.0_idx + (dtdz_up * Ks(i)+ dtdz_dw * Ks(i-1))         ! Coefficient of S(i)_new
       C(i) = - 1.0_idx * dtdz_up * Ks(i)                 ! Coefficient of S(i+1)_new
       D(i) = salt(i) + dt * salt_correct(i)         ! Coefficient of RHS
    end do
    ! Boundary Condition--------------------------------
    ! At the surface---------------
    dtdz_dw = dt / ((zz(N+1)-zz(N))*(depth(N)-depth(N-1)))
    A(N) = -1.0_idx * dtdz_dw * ks(N-1)    ! Coefficient of S(N-1)
    B(N) = 1.0_idx + dtdz_dw * ks(N-1)     ! Coefficient of S(N)
    C(N) = 0.0_idx
    D(N) = salt(N) + dt * ssflux / (zz(N+1)-zz(N)) + dt * salt_correct(N) ! KS(N-1)*dS/dz = S(N) * (E-p) => KS(N-1)*(S(N)-S(N-1)) / dz = (E-P)*S(N)
    ! At the bottom----------------
    A(1) = 0.0_idx     ! Coefficient of S(0)
    B(1) = -1.0_idx*ks(1)/ (zz(2)-zz(1))          ! Coefficient of S(1)
    C(1) = ks(1) / (zz(2)-zz(1))
    D(1) =  0.0_idx            ! KS(1)*dS/dz=0.0
    !---------------------------------------------------
    ! Solve
    call solve_tri_implicit_real_nolap(N,A,B,C,D,Salt_next)
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
  end subroutine cal_next_salt
  !---------------------------------! 
  ! Velocity                        !
  !---------------------------------!
  ! subroutine for calculating new velocity
  subroutine cal_next_uv(dt,N,depth,zz,u,v,km,f,tau_x,tau_y,u_correct,v_correct,u_next,v_next)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: dt,depth(n),zz(n+1),u(n),v(n),km(n-1)
    real(idx),intent(in) :: f,tau_x,tau_y
    real(idx),intent(in) :: u_correct(n),v_correct(n)
    real(idx),intent(inout) :: u_next(n),v_next(n)
    complex(kind(0d0)) :: uv(n),uv_correct(n),A(N),B(N),C(N),D(N),uv_next(N)
    real(idx) :: dtdz_dw,dtdz_up
    integer :: i
    uv=0.0_idx ; uv_next = 0.0_idx; uv_correct=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx
    ! Make matrix
    do i = 2,N-1
       uv(i) = u(i) + (0.0_idx,1.0_idx) * v(i)
       uv_correct(i) = u_correct(i) + (0.0_idx,1.0_idx) * v_correct(i)
       dtdz_dw = dt / ((zz(i+1)-zz(i))*(depth(i)-depth(i-1)))
       dtdz_up = dt / ((zz(i+1)-zz(i))*(depth(i+1)-depth(i)))
       A(i) = cmplx(-1.0_idx * dtdz_dw * Km(i-1),0.0_idx)               ! Coefficient of U(i-1)_new
       B(i) = cmplx(1.0_idx + (dtdz_dw*Km(i-1)+dtdz_up*Km(i)),f*dt)        ! Coefficient of U(i)_new
       C(i) = cmplx(-1.0_idx * dtdz_up * Km(i),0.0_idx)             ! Coefficient of U(i+1)_new
       D(i) = uv(i)+ dt * uv_correct(i)               ! Coefficient of RHS
    end do
    ! Boundary Condition
    ! At the surface
    uv(N) = u(N) + (0.0_idx,1.0_idx) * v(N)
    uv_correct(N) = u_correct(N) + (0.0_idx,1.0_idx) * v_correct(N)

    dtdz_dw = dt / ((zz(N+1)-zz(N))*(depth(N)-depth(N-1)))
    A(N) = cmplx(-1.0_idx * dtdz_dw * Km(N-1),0.0_idx)   ! Coefficient of U,V(N-1)        
    B(N) = cmplx(1.0_idx + dtdz_dw * Km(N-1),f*dt)       ! Coefficient of U,V(N)
    C(N) = cmplx(0.0_idx,0_idx)                   ! Coefficient of U,V(N+1)<-does not exist
    D(N) = uv(N) +  dt * cmplx(tau_x/rho,tau_y/rho) / (zz(N+1)-zz(N))+ dt * uv_correct(N)
    ! At the bottom
    A(1) = cmplx(0.0_idx,0.0_idx)                              ! Coefficient of U,V(0)<-does not exist
    B(1) = cmplx(-1.0_idx*Km(1) / (zz(2)-zz(1)),0.0_idx)       ! Coefficient of U,V(1)
    C(1) = cmplx(Km(1) / (zz(2)-zz(1)),0.0_idx)            ! Coefficient of U,V(2)
    D(1) = cmplx(0.0_idx,0.0_idx)
    call solve_tri_implicit_complex_nolap(N,A,B,C,D,uv_next)
    u_next=real(uv_next,idx) ; v_next=imag(uv_next)
  end subroutine cal_next_uv

  !---------------------------------
  ! TKE
  !---------------------------------
  ! subroutine for calculating new tke
  subroutine cal_next_qq(dt,N,depth,zz,qq,bvf,shear2,L,kq,Km,Kt,tau_x,tau_y,B_1,qq_next)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: dt,depth(N),zz(N+1),qq(N+1),kq(N),km(N-1),kt(N-1)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1),L(N+1),tau_x,tau_y
    real(idx),intent(in) :: B_1
    real(idx),intent(inout) :: qq_next(N+1)
    real(idx),allocatable :: A(:),B(:),C(:),D(:)
    real(idx),allocatable :: RHS_qq(:),SP(:),BP(:),DS(:)
    integer :: i
    real(idx) :: dtdz_dw,dtdz_up,tau_ts
    allocate(A(1:N+1));allocate(B(1:N+1));allocate(C(1:N+1));allocate(D(1:N+1))
    allocate(RHS_qq(1:N+1));allocate(SP(N+1));allocate(BP(N+1));allocate(DS(N+1))
    RHS_qq = 0.0_idx ; SP = 0.0_idx ; BP=0.0_idx
    A=0.0_idx ; B=0.0_idx ; C=0.0_idx ; D=0.0_idx

    ! qq(1)       qq(2)........qq(N)...qq(N+1) 
    !       T(1)                     T(N)
    !             kt(1)        kt(N-1)
    !       kq(1)                    kq(N)
    !             bvf(1)       bvf(n-1)
    ! Make matrix
    do i = 2,N
       dtdz_dw = dt / ((depth(i)-depth(i-1))*(zz(i)-zz(i-1)))
       dtdz_up = dt / ((depth(i)-depth(i-1))*(zz(i+1)-zz(i)))
       A(i) = - 1.0_idx * dtdz_dw * kq(i-1)                     ! Coefficient of qq(i-1)_new
       B(i) = 1.0_idx + (dtdz_dw * kq(i-1) + dtdz_up * kq(i))   ! Coefficient of qq(i)_new
       C(i) = - 1.0_idx * (dtdz_up * kq(i))                ! Coefficient of qq(i+1)_new
       !------------------------------------------------
       ! Compute Shear production and Buoyancy production
       !------------------------------------------------
       ! Buoyancy frequency
       SP(i) = km(i-1) * Shear2(i-1)
       BP(i) = -1.0_idx * kt(i-1) * bvf(i-1)
       DS(i) = (qq(i)**1.5_idx) / (B_1 * l(i))
       RHS_qq(i) = dt * 2.0_idx * (SP(i)+BP(i)-DS(i)) ! Forcing term
       D(i) = qq(i) +  RHS_qq(i)                  ! Coefficient of RHS
    end do
    ! Boundary Condition---------------------------------
    ! At the surface-----------
    tau_ts=sqrt(tau_x**2+tau_y**2) / rho ! = u_star^2 (tau / rho)
    A(N+1) = 0.0_idx ! Coefficient of qq_new(N) 
    B(N+1) = 1.0_idx ! Coefficient of qq_new(N+1)
    C(N+1) = 0.0_idx ! Coefficient of qq_new(N+2) <-does not exist
    D(N+1) = ((B_1)**(2.0_idx/3.0_idx))*(tau_ts)! surface energy
    ! At the bottom------------
    A(1) = 0.0_idx   ! Coefficient of qq(0)_new
    B(1) = 1.0_idx   ! Coefficient of qq(1)_new
    C(1) = 0.0_idx   ! Coefficient of qq(2)_new
    D(1) = 0.0_idx
    !----------------------------------------------------
    ! Solve---------------------------------------------
    call solve_tri_implicit_real_nolap(N+1,A,B,C,D,qq_next)
    do i=1,N+1
!       write(*,*) qq_next(i)
       if(qq_next(i) .le. qq_min) then
          qq_next(i) = qq_min
       end if
    end do
    deallocate(A);deallocate(B);deallocate(C);deallocate(D)
    deallocate(RHS_qq);deallocate(SP);deallocate(BP);deallocate(DS)
  end subroutine cal_next_qq
end module solve_diag
