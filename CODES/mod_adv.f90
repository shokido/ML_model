module mod_adv
  !================================================
  ! Module program for external advection file
  !================================================
  use param
  implicit none
  !--------------------------------
  ! Temperature
  logical :: use_adv_temp
  character(maxlen),allocatable :: fnames_adv_temp(:)
  character(maxlen) :: varname_adv_temp
  integer :: ntime_adv_temp,nfile_adv_temp
  integer :: ind1_adv_temp,ind2_adv_temp
  real(idx) :: wgt1_adv_temp,wgt2_adv_temp
  real(idx),allocatable :: time_adv_temp(:)
  real(idx),allocatable :: temp_adv(:,:,:,:)
  real(idx),allocatable :: temp_adv_1d(:)
  ! Salinity
  logical :: use_adv_salt
  character(maxlen),allocatable :: fnames_adv_salt(:)
  character(maxlen) :: varname_adv_salt
  integer :: ntime_adv_salt,nfile_adv_salt
  integer :: ind1_adv_salt,ind2_adv_salt
  real(idx) :: wgt1_adv_salt,wgt2_adv_salt
  real(idx),allocatable :: time_adv_salt(:)
  real(idx),allocatable :: salt_adv(:,:,:,:)
  real(idx),allocatable :: salt_adv_1d(:)
  ! U
  logical :: use_adv_u
  character(maxlen),allocatable :: fnames_adv_u(:)
  character(maxlen) :: varname_adv_u
  integer :: ntime_adv_u,nfile_adv_u
  integer :: ind1_adv_u,ind2_adv_u
  real(idx) :: wgt1_adv_u,wgt2_adv_u
  real(idx),allocatable :: time_adv_u(:)
  real(idx),allocatable :: u_adv(:,:,:,:)
  real(idx),allocatable :: u_adv_1d(:)
  ! V
  logical :: use_adv_v
  character(maxlen),allocatable :: fnames_adv_v(:)
  character(maxlen) ::  varname_adv_v
  integer :: ntime_adv_v,nfile_adv_v
  integer :: ind1_adv_v,ind2_adv_v
  real(idx) :: wgt1_adv_v,wgt2_adv_v
  real(idx),allocatable :: time_adv_v(:)
  real(idx),allocatable :: v_adv(:,:,:,:)
  real(idx),allocatable :: v_adv_1d(:)
  ! 
  namelist/input_adv_param/use_adv_temp,nfile_adv_temp,varname_adv_temp
  namelist/input_adv_param/use_adv_salt,nfile_adv_salt,varname_adv_salt
  namelist/input_adv_param/use_adv_u,nfile_adv_u,varname_adv_u
  namelist/input_adv_param/use_adv_v,nfile_adv_v,varname_adv_v
  namelist/input_adv_io/fnames_adv_temp
  namelist/input_adv_io/fnames_adv_salt
  namelist/input_adv_io/fnames_adv_u
  namelist/input_adv_io/fnames_adv_v
end module mod_adv
