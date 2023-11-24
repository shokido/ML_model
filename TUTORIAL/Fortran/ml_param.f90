module ml_param
  implicit none
  ! Select float (idx=4) or double (idx=8)
  integer,parameter :: idx=8
  ! Maximum length of character
  integer,parameter :: maxlen=400
  real(idx),parameter :: missing_value=9999.0_idx
  real(idx),parameter :: rho=1.024e3_idx   ! reference density of seawater
  real(idx),parameter :: cpw=3986.0_idx,g=9.8_idx,T0=20.0_idx,S0=34.0_idx
  real(idx),parameter :: qq_min=1.0e-8
  real(idx),parameter :: day_to_sec=60.0_idx*60.0_idx*24.0_idx
  ! Parameter setting
  real(idx) :: nu                   ! background diffusion
  real(idx) :: nu_t                 ! backgroud diffusion of temperature
  real(idx) :: nu_s                 ! backgroud diffusion of salinity
  real(idx) :: beta1,beta2,r_long   ! shortwave penetration parameters
end module ml_param
