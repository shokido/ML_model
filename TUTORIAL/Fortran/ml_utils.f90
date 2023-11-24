module ml_utils
  use ml_param
  implicit none
  private
  type datetime
     integer :: year,month,day,hour,minute,second
  end type datetime
  real(idx),parameter :: kappa=0.41_idx  
  real(idx),parameter :: q2_init = 1.0e-8_idx
  real(idx) :: tiny=1.0e-12_idx
  public :: cal_pdens,cal_alpha,cal_beta,cal_bvf,cal_shear
  public :: heat_absorb,cal_penet
  public :: initialize_turb
  public :: kappa,cal_mo_inv
  public :: datetime,dt_to_jd,jd_to_dt
contains
  !============================================================================
  function cal_pdens(s,t) result(ret)
    implicit none
    real(idx),intent(in) :: s,t
    real(idx) :: rho_w,ret
    rho_w = 999.842594_idx + t * 6.793952e-2_idx - (t**2) *  9.095290e-3_idx &
         &          + 1.001685e-04_idx * t**3 - 1.120083e-6_idx * t**4 + 6.536332e-9_idx * t**5 
    ret = rho_w  + (0.824493_idx-4.0899e-3_idx * t + (7.6438e-5_idx) * t**2 -8.2467e-7_idx &
         &               * t**3 + 5.3875e-9_idx * t**4)*s &
         &                + (-5.72466e-3_idx + 1.0227e-4_idx * t - 1.6546e-6_idx * t**2) * s**(1.5_idx) &
         &                + 4.8314e-4_idx * s**2
  end function cal_pdens
  !============================================================================
  !-Calculate thermal expansion coefficient
  ! Note that this function generally returns negative value
  function cal_alpha(s,t) result(alpha)
    implicit none
    real(idx),intent(in) :: s,t
    real(idx) :: rho_1,rho_2,dt
    real(idx) :: alpha
    dt=0.01_idx
    rho_1=cal_pdens(s,t)
    rho_2=cal_pdens(s,t+dt)
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
    rho_1=cal_pdens(s+ds,t)
    rho_2=cal_pdens(s,t)    ! rho_1 > rho_2
    beta = (rho_1-rho_2)/(rho_1*ds)  ! beta > 0
  end function cal_beta
  !============================================================================
  function cal_bvf(nz,z_rho,temp,salt) result(bv)
    implicit none
    integer,intent(in) :: nz
    real(idx),intent(in) :: z_rho(nz),temp(nz),salt(nz)
    real(idx) :: bv(nz-1)
    integer :: iz
    do iz=1,nz-1
       ! Buoyancy frequency
       bv(iz) =  - 1.0_idx * g * (cal_pdens(salt(iz),temp(iz))- cal_pdens(salt(iz+1),temp(iz+1))) / &
            & ((z_rho(iz)-z_rho(iz+1))* &
            & 0.5*(cal_pdens(salt(iz),temp(iz))+cal_pdens(salt(iz+1),temp(iz+1))))
    end do
  end function cal_bvf
  function cal_shear(nz,z_rho,u,v) result(shear)
    implicit none
    integer,intent(in) :: nz
    real(idx),intent(in) :: z_rho(nz),u(nz),v(nz)
    real(idx) :: shear(nz-1)
    real(idx) :: uz,vz
    integer :: iz
    do iz=1,nz-1
       ! Shear
       uz = (u(iz)-u(iz+1)) / (z_rho(iz)-z_rho(iz+1))
       vz = (v(iz)-v(iz+1)) / (z_rho(iz)-z_rho(iz+1))
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
  function cal_mo_inv(sst,sss,tau_x,tau_y,sw,nsw,salflux) result(mo_inv)
    real(idx),intent(in) :: sst,sss,tau_x,tau_y,sw,nsw,salflux
    real(idx) :: mo_inv
    real(idx) :: alpha,beta,bf,u_star2,u_star3
    alpha=cal_alpha(sss,sst) ! < 0
    beta=cal_beta(sss,sst)   ! > 0
    !bf=-wb0=-g*alpha*wto+beta*ws0
    bf = -1.0_idx * g * (alpha * (sw+nsw) / (rho*cpw)+beta*salflux)
    ! tau_x = rho * u_star^2 * (direction)
    ! tau_y = rho * v_star^2
    ! u_star3 = (sqrt(u_star^2+v_star^2))^3
    u_star2=sqrt(tau_x**2+tau_y**2)/rho
    u_star3 = sqrt(u_star2)**3
    !    mo =u_star3 / (dkappa * bf)
    u_star3=max(u_star3,1.0e-10)
    mo_inv = kappa * bf / u_star3
    !bf <0 -> mo < 0 ! unstable
    !bf >0 -> mo > 0 ! stable
  end function cal_mo_inv
  function dt_to_jd(dt) result(jd)
    implicit none
    type(datetime),intent(in) :: dt
    integer :: year,month,day
    integer :: a,y,m
    real(idx) :: jd
    year = dt%year;month=dt%month;day=dt%day
    a = int((14-month)/12)
    y = year + 4800 - a
    m = month + 12 * a - 3
    ! 12:00's JDN
    jd = day + int((153*m+2)/5) + 365 * y + int(y/4) - int(y/100) &
         & + int(y/400)-32045+dt%hour/24.0_idx+dt%minute/(24.0_idx*60.0_idx)&
         & +dt%second/(24.0_idx*60.0_idx*60.0_idx)
  end function dt_to_jd
  function jd_to_dt(jd) result(dt)
    implicit none
    real(idx),intent(in) :: jd
    integer:: year,month,day
    integer:: hour,minute,second
    integer :: f,e,g,h
    real(idx) :: res
    type(datetime) :: dt
    res=jd-int(jd)
    f = int(jd) + 1401 + (((4 * int(jd) + 274277) / 146097) * 3) / 4 - 38
    e = 4 * f + 3
    g = (mod(e,1461)/4)
    h = 5 * g + 2
    day = (mod(h,153)/5)+1
    month = mod((h/153)+2,12)+1
    year = (e/1461) - 4716 + (12 + 2 - month) / 12
    dt%year=year;dt%month=month;dt%day=day
    hour=int(res*24)
    res=res-hour/(24.0_idx)
    minute=int(res*60*24)
    res=res-minute/(60.0_idx*24.0_idx)
    second=int(res*60*24*60)
    dt%hour=hour;dt%minute=minute;dt%second=second
  end function jd_to_dt
end module ml_utils


