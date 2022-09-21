module mod_clm
  !================================================
  ! Module program for external climatology file
  !================================================
  use param
  implicit none
  !--------------------------------
  ! Temperature
  logical :: use_clm_temp
  character(maxlen),allocatable :: fnames_clm_temp(:)
  character(maxlen) :: varname_clm_temp
  integer :: ntime_clm_temp,nfile_clm_temp
  integer :: ind1_clm_temp,ind2_clm_temp
  real(idx) :: wgt1_clm_temp,wgt2_clm_temp
  real(idx),allocatable :: time_clm_temp(:)
  real(idx),allocatable :: temp_clm(:,:,:,:)
  real(idx),allocatable :: temp_clm_1d(:)
  ! Salinity
  logical :: use_clm_salt
  character(maxlen),allocatable :: fnames_clm_salt(:)
  character(maxlen) :: varname_clm_salt
  integer :: ntime_clm_salt,nfile_clm_salt
  integer :: ind1_clm_salt,ind2_clm_salt
  real(idx) :: wgt1_clm_salt,wgt2_clm_salt
  real(idx),allocatable :: time_clm_salt(:)
  real(idx),allocatable :: salt_clm(:,:,:,:)
  real(idx),allocatable :: salt_clm_1d(:)
  ! U
  logical :: use_clm_u
  character(maxlen),allocatable :: fnames_clm_u(:)
  character(maxlen) :: varname_clm_u
  integer :: ntime_clm_u,nfile_clm_u
  integer :: ind1_clm_u,ind2_clm_u
  real(idx) :: wgt1_clm_u,wgt2_clm_u
  real(idx),allocatable :: time_clm_u(:)
  real(idx),allocatable :: u_clm(:,:,:,:)
  real(idx),allocatable :: u_clm_1d(:)
  ! V
  logical :: use_clm_v
  character(maxlen),allocatable :: fnames_clm_v(:)
  character(maxlen) ::  varname_clm_v
  integer :: ntime_clm_v,nfile_clm_v
  integer :: ind1_clm_v,ind2_clm_v
  real(idx) :: wgt1_clm_v,wgt2_clm_v
  real(idx),allocatable :: time_clm_v(:)
  real(idx),allocatable :: v_clm(:,:,:,:)
  real(idx),allocatable :: v_clm_1d(:)
  ! 
  namelist/input_clm_param/use_clm_temp,nfile_clm_temp,varname_clm_temp
  namelist/input_clm_param/use_clm_salt,nfile_clm_salt,varname_clm_salt
  namelist/input_clm_param/use_clm_u,nfile_clm_u,varname_clm_u
  namelist/input_clm_param/use_clm_v,nfile_clm_v,varname_clm_v
  namelist/input_clm_io/fnames_clm_temp
  namelist/input_clm_io/fnames_clm_salt
  namelist/input_clm_io/fnames_clm_u
  namelist/input_clm_io/fnames_clm_v
end module mod_clm
