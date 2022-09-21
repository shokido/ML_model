module mod_rst
  !================================================
  ! Module program for create/write restart file
  !================================================
  use param
  use calendar_sub
  use ncdf_write
  implicit none
  ! 3D Oceanic array
  !  integer :: out_rst_flag,out_rst_int
  character(maxlen) :: fname_out_rst_rho
  real(idx),allocatable :: temp_rst(:,:,:),salt_rst(:,:,:),u_rst(:,:,:),v_rst(:,:,:)
  character(maxlen) :: fname_out_rst_q
  real(idx),allocatable :: qq_rst(:,:,:),l_rst(:,:,:)
  ! Namelist
  namelist/output_rst_flag/fname_out_rst_rho
  namelist/output_rst_flag/fname_out_rst_q
  !========================================================
contains
  subroutine prepare_output_rst_rho(fname,nlon,nlat,nlev,lon,lat,lev,time_in, &
       & ref_yymmdd,ref_hhmmss,flag_yymmdd,missing_value)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    integer,parameter :: ntime=1
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev),time_in
    real(idx) :: time(ntime)
    integer,intent(in) :: ref_yymmdd,ref_hhmmss,flag_yymmdd
    real(idx),intent(in) :: missing_value
    character :: time_unit*100
    time_unit=calendar_create_time_att(ref_yymmdd,ref_hhmmss,flag_yymmdd)
    time(1)=time_in
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","z_rho","time",&
         & "degrees_east","degrees_north","m",trim(time_unit), &
         & lon,lat,lev,time)
    call writenet_dv(trim(fname),"lon","lat","z_rho","time",&
         & 1,(/"temp"/),(/"degrees celcius"/),missing_value)
    call writenet_dv(trim(fname),"lon","lat","z_rho","time",&
         & 1,(/"salt"/),(/"psu"/),missing_value)
    call writenet_dv(trim(fname),"lon","lat","z_rho","time",&
         & 1,(/"u"/),(/"m/s"/),missing_value)
    call writenet_dv(trim(fname),"lon","lat","z_rho","time",&
         & 1,(/"v"/),(/"m/s"/),missing_value)
  end subroutine prepare_output_rst_rho
  ! Output q-lev file
  subroutine prepare_output_rst_q(fname,nlon,nlat,nlev,lon,lat,z_q,time_in, &
       & ref_yymmdd,ref_hhmmss,flag_yymmdd,missing_value)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: lon(nlon),lat(nlat),z_q(nlev),time_in
    integer,intent(in) :: ref_yymmdd,ref_hhmmss,flag_yymmdd
    integer,parameter :: ntime=1
    real(idx) :: time(ntime)
    real(idx) :: missing_value
    character :: time_unit*100
    time_unit=calendar_create_time_att(ref_yymmdd,ref_hhmmss,flag_yymmdd)
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","z_q","time",&
         & "degrees_east","degrees_north","m",trim(time_unit), &
         & lon,lat,z_q,time)
    call writenet_dv(trim(fname),"lon","lat","z_q","time",&
         & 1,(/"qq"/),(/"m^2/s^2"/),missing_value)
    call writenet_dv(trim(fname),"lon","lat","z_q","time",&
            & 1,(/"l"/),(/"m"/),missing_value)
  end subroutine prepare_output_rst_q
  subroutine output_rst_rho(fname,nlon,nlat,nlev,temp,salt,u,v)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: temp(nlon,nlat,nlev),salt(nlon,nlat,nlev),u(nlon,nlat,nlev),v(nlon,nlat,nlev)
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    tmp_4D=0.0_idx
    tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp(1:nlon,1:nlat,1:nlev)
    call writenet_wv(trim(fname),"temp",(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp_4D)
    tmp_4D(1:nlon,1:nlat,1:nlev,1)=salt(1:nlon,1:nlat,1:nlev)
    call writenet_wv(trim(fname),"salt",(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp_4D)
    tmp_4D(1:nlon,1:nlat,1:nlev,1)=u(1:nlon,1:nlat,1:nlev)
    call writenet_wv(trim(fname),"u",(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp_4D)
    tmp_4D(1:nlon,1:nlat,1:nlev,1)=v(1:nlon,1:nlat,1:nlev)
    call writenet_wv(trim(fname),"v",(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp_4D)
  end subroutine output_rst_rho
  subroutine output_rst_q(fname,nlon,nlat,nlev,qq,l)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: qq(nlon,nlat,nlev),l(nlon,nlat,nlev)
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    tmp_4D=0.0_idx
    tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq(1:nlon,1:nlat,1:nlev)
    call writenet_wv(trim(fname),"qq",(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp_4D)
    tmp_4D(1:nlon,1:nlat,1:nlev,1)=l(1:nlon,1:nlat,1:nlev)
    call writenet_wv(trim(fname),"l",(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp_4D)
  end subroutine output_rst_q
end module mod_rst
