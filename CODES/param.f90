module param
  implicit none
  ! Select float (idx=4) or double (idx=8)
  integer,parameter :: idx=8
  ! Maximum length of character
  integer,parameter :: maxlen=400
  real(idx),parameter :: missing_value=9999.0_idx
  ! Namelist---------------------------------------------------------
  integer,parameter :: nmlf_master=10,nmlf_set=20,nmlf_io=30
  character(maxlen) :: namelist_master
  character(maxlen) :: namelist_set,namelist_io
  ! Time step========================================================
  real(idx) :: dt
  integer :: start_yymmdd_int,start_hhmmss_int
  integer :: end_yymmdd_int,end_hhmmss_int
  real(idx),parameter :: sec_to_day = 1.0_idx / (24.0*60.0*60.0_idx)
  real(idx),parameter :: day_to_sec = 60.0*60.0*24.0_idx
  real(idx),parameter :: year_to_sec=60.0_idx * 60.0_idx * 24.0_idx * 365 ! [s/year]
  ! Parameters
  real(idx),parameter :: rho=1.024e3_idx   ! reference density of seawater
  real(idx),parameter :: cpw=3986.0_idx,g=9.8_idx,T0=20.0_idx,S0=34.0_idx
  !==================================================================
  ! Parameter setting
  real(idx) :: nu                   ! background diffusion
  real(idx) :: nu_t                 ! backgroud diffusion of temperature
  real(idx) :: nu_s                 ! backgroud diffusion of salinity
  real(idx) :: beta1,beta2,r_long   ! shortwave penetration parameters
  real(idx) :: tau_temp,tau_salt,tau_u,tau_v
  real(idx) :: tau_temp_day,tau_salt_day,tau_u_day,tau_v_day
  namelist/master/namelist_set,namelist_io
  namelist/time/dt,start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int
  namelist/parameters/nu,nu_t,nu_s
  namelist/parameters/beta1,beta2,r_long
  namelist/parameters/tau_temp_day,tau_salt_day,tau_u_day,tau_v_day
end module param
