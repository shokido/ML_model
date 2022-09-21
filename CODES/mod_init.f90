module mod_init
  !================================================
  ! Module program for read initital condition
  !================================================
  use param
  implicit none
  !--------------------------------
  ! Temperature
  character(maxlen) :: fname_init_temp,varname_init_temp
  integer :: istart_temp
  ! Salinity
  character(maxlen) :: fname_init_salt,varname_init_salt
  integer :: istart_salt
  ! U
  character(maxlen) :: fname_init_u,varname_init_u
  integer :: istart_u
  ! V
  character(maxlen) :: fname_init_v,varname_init_v
  integer :: istart_v
  ! QQ
  character(maxlen) :: fname_init_qq,varname_init_qq
  integer :: istart_qq
  ! L
  character(maxlen) :: fname_init_l,varname_init_l
  integer :: istart_l
  
  namelist/input_init/fname_init_temp,varname_init_temp,istart_temp
  namelist/input_init/fname_init_salt,varname_init_salt,istart_salt
  namelist/input_init/fname_init_u,varname_init_u,istart_u
  namelist/input_init/fname_init_v,varname_init_v,istart_v
  namelist/input_init/fname_init_qq,varname_init_qq,istart_qq
  namelist/input_init/fname_init_l,varname_init_l,istart_l
end module mod_init
