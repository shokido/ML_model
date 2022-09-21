module mod_lin
  !================================================
  ! Module program for linear operator
  !================================================
  use param
  use ncdf_write
  implicit none
  character(maxlen) :: fname_out_lin_rho
  real(idx),allocatable :: matrix(:,:,:,:)
  ! Namelist
  namelist/output_lin_flag/fname_out_lin_rho
  !========================================================
contains
  subroutine prepare_output_lin(fname,nlon,nlat,nlev,lon,lat,lev,missing_value)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    integer,parameter :: ntime=1
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev)
    real(idx),intent(in) :: missing_value
    real(idx) :: lev_new(nlev*4)
    lev_new(1:nlev)=lev(1:nlev)
    lev_new(nlev+1:nlev*2)=lev(1:nlev)*2
    lev_new(2*nlev+1:nlev*3)=lev(1:nlev)*3
    lev_new(3*nlev+1:nlev*4)=lev(1:nlev)*4
    call writenet_pre(fname,nlon,nlat,nlev*4,nlev*4,"lon","lat","z_rho1","z_rho2",&
         & "degrees_east","degrees_north","m","m", &
         & lon,lat,lev_new,lev_new)
    call writenet_dv(trim(fname),"lon","lat","z_rho1","z_rho2",&
         & 1,(/"matrix"/),(/""/),missing_value)
  end subroutine prepare_output_lin
  subroutine output_lin_rho(fname,nlon,nlat,nlev,mat)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: mat(nlon,nlat,nlev*4,nlev*4)
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:4*nlev,1:4*nlev)
    tmp_4D=0.0_idx
    tmp_4D(1:nlon,1:nlat,1:4*nlev,1:4*nlev)=mat(1:nlon,1:nlat,1:4*nlev,1:4*nlev)
    call writenet_wv(trim(fname),"matrix",1,nlon,1,nlat,1,4*nlev,1,4*nlev,tmp_4D)
  end subroutine output_lin_rho
end module mod_lin
