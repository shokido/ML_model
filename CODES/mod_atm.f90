module mod_atm
#include 'CPPLISTS.h'
  !================================================
  ! Module program for external atmospheric
  !================================================
  use param
  implicit none
  !--------------------------------
  ! Precipitation
  character(maxlen),allocatable :: fnames_pr(:)
  character(maxlen) :: varname_pr
  integer :: ntime_pr,nfile_pr
  integer :: ind1_pr,ind2_pr
  real(idx) :: wgt1_pr,wgt2_pr
  real(idx),allocatable :: time_pr(:)
  ! SW radiation
  character(maxlen),allocatable :: fnames_sw(:)
  character(maxlen) :: varname_sw
  integer :: ntime_sw,nfile_sw
  integer :: ind1_sw,ind2_sw
  real(idx) :: wgt1_sw,wgt2_sw
  real(idx),allocatable :: time_sw(:)
  ! LW radiation
  character(maxlen),allocatable :: fnames_lw(:)
  character(maxlen) :: varname_lw
  integer :: ntime_lw,nfile_lw
  integer :: ind1_lw,ind2_lw
  real(idx) :: wgt1_lw,wgt2_lw
  real(idx),allocatable :: time_lw(:)
  ! Surface air/Sensible heat
  character(maxlen),allocatable :: fnames_ta(:)
  character(maxlen) :: varname_ta
  integer :: ntime_ta,nfile_ta
  integer :: ind1_ta,ind2_ta
  real(idx) :: wgt1_ta,wgt2_ta
  real(idx),allocatable :: time_ta(:)
  ! Humidity/Latent heat
  character(maxlen),allocatable :: fnames_qa(:)
  character(maxlen) :: varname_qa
  integer :: ntime_qa,nfile_qa
  integer :: ind1_qa,ind2_qa
  real(idx) :: wgt1_qa,wgt2_qa
  real(idx),allocatable :: time_qa(:)
  ! Zonal wind speed/ wind stress
  character(maxlen),allocatable :: fnames_uw(:)
  character(maxlen) :: varname_uw
  integer :: ntime_uw,nfile_uw
  integer :: ind1_uw,ind2_uw
  real(idx) :: wgt1_uw,wgt2_uw
  real(idx),allocatable :: time_uw(:)
  ! Meridional wind speed/wind stress
  character(maxlen),allocatable :: fnames_vw(:)
  character(maxlen) :: varname_vw
  integer :: ntime_vw,nfile_vw
  integer :: ind1_vw,ind2_vw
  real(idx) :: wgt1_vw,wgt2_vw
  real(idx),allocatable :: time_vw(:)
#if defined(WSPEED_IN)
  ! Wind speed
  character(maxlen),allocatable :: fnames_ws(:)
  character(maxlen) :: varname_ws
  integer :: ntime_ws,nfile_ws
  integer :: ind1_ws,ind2_ws
  real(idx) :: wgt1_ws,wgt2_ws
  real(idx),allocatable :: time_ws(:)
#endif

  real(idx) :: ta,qa
  real(idx) :: uwind,vwind
  real(idx) :: ws
  real(idx) :: ld
  real(idx) :: heat_solar,heat_no_solar,e_p,ssflux
  real(idx),allocatable :: sw_in(:,:,:),lw_in(:,:,:),pr_in(:,:,:)
  real(idx),allocatable :: ta_in(:,:,:),qa_in(:,:,:)
  real(idx),allocatable :: uw_in(:,:,:),vw_in(:,:,:)
#if defined(WSPEED_IN)
  real(idx),allocatable :: ws_in(:,:,:)
#endif
  !namelist/input_atm/fnames_atm
  namelist/input_atm_param/nfile_pr,varname_pr
  namelist/input_atm_param/nfile_sw,varname_sw
  namelist/input_atm_param/nfile_lw,varname_lw
  namelist/input_atm_param/nfile_ta,varname_ta
  namelist/input_atm_param/nfile_qa,varname_qa
  namelist/input_atm_param/nfile_uw,varname_uw
  namelist/input_atm_param/nfile_vw,varname_vw
#if defined(WSPEED_IN)
  namelist/input_atm_param/nfile_ws,varname_ws
#endif

  namelist/input_atm_io/fnames_pr
  namelist/input_atm_io/fnames_sw
  namelist/input_atm_io/fnames_lw
  namelist/input_atm_io/fnames_ta
  namelist/input_atm_io/fnames_qa
  namelist/input_atm_io/fnames_uw
  namelist/input_atm_io/fnames_vw
#if defined(WSPEED_IN)
  namelist/input_atm_io/fnames_ws
#endif
end module mod_atm
