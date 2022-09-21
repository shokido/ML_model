module eos_sub
  use param
  implicit none
  private
  real(idx),parameter :: kappa=0.41_idx  
  real(idx),parameter :: q2_init = 1.0e-8_idx
  real(idx) :: tiny=1.0e-12_idx
  public :: cal_dens,cal_alpha,cal_beta,cal_bv,cal_shear
  public :: heat_absorb,cal_penet
  public :: initialize_turb
  public :: kappa
contains
  !============================================================================
  function cal_dens(s,t) result(ret)
    implicit none
    real(idx),intent(in) :: s,t
    real(idx) :: rho_w,ret
    rho_w = 999.842594_idx + t * 6.793952e-2_idx - (t**2) *  9.095290e-3_idx &
         &          + 1.001685e-04_idx * t**3 - 1.120083e-6_idx * t**4 + 6.536332e-9_idx * t**5 
    ret = rho_w  + (0.824493_idx-4.0899e-3_idx * t + (7.6438e-5_idx) * t**2 -8.2467e-7_idx &
         &               * t**3 + 5.3875e-9_idx * t**4)*s &
         &                + (-5.72466e-3_idx + 1.0227e-4_idx * t - 1.6546e-6_idx * t**2) * s**(1.5_idx) &
         &                + 4.8314e-4_idx * s**2
  end function cal_dens
  !============================================================================
  !-Calculate thermal expansion coefficient
  ! Note that this function generally returns negative value
  function cal_alpha(s,t) result(alpha)
    implicit none
    real(idx),intent(in) :: s,t
    real(idx) :: rho_1,rho_2,dt
    real(idx) :: alpha
    dt=0.01_idx
    rho_1=cal_dens(s,t)
    rho_2=cal_dens(s,t+dt)
    alpha = (rho_2-rho_1)/(rho_1*dt)  ! alpha < 0
  end function cal_alpha
  !============================================================================
  !-Calculate salinity extraction coefficient
  function cal_beta(s,t) result(beta)
    implicit none
    real(idx),intent(in) :: s,t
    real(idx) :: rho_1,rho_2,ds
    real(idx) :: beta
    ds=0.01_idx
    rho_1=cal_dens(s+ds,t)
    rho_2=cal_dens(s,t)    ! rho_1 > rho_2
    beta = (rho_1-rho_2)/(rho_1*ds)  ! beta > 0
  end function cal_beta
  !============================================================================
  function cal_bv(nz,z_rho,temp,salt) result(bv)
    implicit none
    integer,intent(in) :: nz
    real(idx),intent(in) :: z_rho(nz),temp(nz),salt(nz)
    real(idx) :: bv(nz-1)
    integer :: iz
    do iz=1,nz-1
       ! Buoyancy frequency
       bv(iz) =  - 1.0_idx * g * (cal_dens(salt(iz+1),temp(iz+1))- cal_dens(salt(iz),temp(iz))) / &
            & ((z_rho(iz+1)-z_rho(iz))* cal_dens(salt(iz),temp(iz)))
    end do
  end function cal_bv
  function cal_shear(nz,z_rho,u,v) result(shear)
    implicit none
    integer,intent(in) :: nz
    real(idx),intent(in) :: z_rho(nz),u(nz),v(nz)
    real(idx) :: shear(nz-1)
    real(idx) :: uz,vz
    integer :: iz
    do iz=1,nz-1
       ! Shear
       uz = (u(iz+1)-u(iz)) / (z_rho(iz+1)-z_rho(iz))
       vz = (v(iz+1)-v(iz)) / (z_rho(iz+1)-z_rho(iz))
       shear(iz) = (uz**2 + vz**2 + tiny)
    end do
  end function cal_shear
  !================================================
  subroutine initialize_turb(N,z,q2,l)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: z(N)
    real(idx),intent(inout) :: q2(N),l(N)
    real(idx),parameter :: H=20.0_idx
    integer :: i
    ! initial condition=========================================
    do i=1,N
       q2(i) = q2_init*exp(z(i)/H)
       l(i) = kappa * abs(z(i))
    end do
  end subroutine initialize_turb
  !============================================================================
  function cal_penet(level,sw,beta1,beta2,r_long) result(ret)
    implicit none
    real(idx) :: level,sw,beta1,beta2,r_long
    real(idx) :: ret
    real(idx) :: z1,z1b1,z1b2
    z1=abs(level)
    z1b1 = min(z1/beta1,100.0_idx) ; z1b2 = min(z1/beta2,100.0_idx)
    ret = sw * (r_long * exp(-1.0_idx * z1b1) +&
         & (1.0_idx - r_long) * exp(-1.0_idx * z1b2))
  end function cal_penet
  subroutine heat_absorb(nz,level,heat_sl,beta1,beta2,r_long,heat)
    implicit none
    integer,intent(in) :: nz
    real(idx),intent(in) :: level(nz)
    real(idx),intent(in) :: heat_sl
    real(idx),intent(in) :: beta1,beta2,r_long
    real(idx),intent(inout) :: heat(nz)
    integer :: i
    do i=1,nz
       heat(i) = cal_penet(level(i),heat_sl,beta1,beta2,r_long)
    end do
  end subroutine heat_absorb
end module eos_sub


