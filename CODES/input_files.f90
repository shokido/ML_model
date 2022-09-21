module input_files
  use param
  implicit none
  private
  public :: read_grid_file
  public :: read_atmos_file,read_atmos_files
  public :: read_clm_file
  public :: read_adv_file
  public :: read_init_file
  public :: linear_int_one
  public :: time_wgt,set_data
  public :: time_wgt2
  public :: read_3D_files,read_4D_files
contains
  !==================================================
  ! Check subroutine
  !==================================================  
  subroutine check_r(status)
    use netcdf
    implicit none
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
       write(*,*) trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check_r
  !==================================================
  ! Subroutine for get dimension
  !================================================== 
  subroutine get_dimension(ncid, name, dim, dims)
    use netcdf
    implicit none
    integer,           intent(in) :: ncid
    character(len=*),  intent(in) :: name
    integer,           intent(inout) :: dim
    real(idx), allocatable, intent(inout) :: dims(:)
    integer :: err
    integer :: varid, dimid
    call check_r(nf90_inq_dimid(ncid, name, dimid) )
    call check_r(nf90_inquire_dimension(ncid, dimid, len=dim) )
    allocate(dims(dim), stat=err)
    if (err /= 0) print *, name, ": Allocation request denied"
    call check_r( nf90_inq_varid(ncid, name, varid) )
    call check_r( nf90_get_var(ncid, varid, dims) )
  end subroutine get_dimension
  subroutine get_var_units(fname_in,varname,units)
    use netcdf
    implicit none
    character(len=*),  intent(in) :: fname_in,varname
    character(len=*), intent(out) :: units
    integer :: ncid,varid
    call check_r(nf90_open(trim(fname_in),nf90_nowrite,ncid))
    call check_r(nf90_inq_varid(ncid,varname, varid))
    call check_r(nf90_get_att(ncid, varid, 'units', units) )
  end subroutine get_var_units
  !==================================================
  ! Read grid file
  !==================================================  
  subroutine read_grid_file(fname,nlon,nlat,nlev,lon_grd,lat_grd,f,z_rho,z_q)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(inout) :: nlon,nlat,nlev
    real(idx), allocatable, intent(inout) :: lon_grd(:),lat_grd(:),f(:),z_rho(:),z_q(:)
    real(idx),allocatable :: lev_grd(:)
    integer :: ncid,varid,start(4),count(4)
    !Open file
    call check_r( nf90_open(trim(fname),nf90_nowrite,ncid))
    if (allocated(lon_grd) .eqv. .true.) then
       deallocate(lon_grd)
    end if
    if (allocated(lat_grd) .eqv. .true.) then
       deallocate(lat_grd)
    end if
    if (allocated(lev_grd) .eqv. .true.) then
       deallocate(lev_grd)
    end if
    if (allocated(f) .eqv. .true.) then
       deallocate(f)
    end if
    if (allocated(z_rho) .eqv. .true.) then
       deallocate(z_rho)
    end if
    if (allocated(z_q) .eqv. .true.) then
       deallocate(z_q)
    end if
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lev",nlev,lev_grd)
    ! Obtain f
    allocate(f(nlat))
    call check_r(nf90_inq_varid(ncid,"f", varid) )
    call check_r(nf90_get_var(ncid, varid,f))
    ! Obtain z_rho
    allocate(z_rho(nlev))
    call check_r( nf90_inq_varid(ncid,"z_rho", varid) )
    call check_r( nf90_get_var(ncid, varid,z_rho) )
    ! Obtain z_q
    allocate(z_q(nlev))
    call check_r( nf90_inq_varid(ncid,"z_q", varid) )
    call check_r( nf90_get_var(ncid, varid,z_q) )
  end subroutine read_grid_file
  function modify_time(ntime,time_in,time_units,start_yymmdd,start_hhmmss) result(time)
    use calendar_sub
    implicit none
    integer,intent(in) :: ntime
    real(idx),intent(in) :: time_in(ntime)
    character(len=*),intent(in) :: time_units
    integer,intent(in) :: start_yymmdd,start_hhmmss
    real(idx) :: time(ntime)
    character :: flag_char*8,yr_char*4,mn_char*2,dy_char*2,hr_char*2,min_char*2,sec_char
    integer :: ref_year,ref_month,ref_day,ref_hour,ref_min,ref_sec
    integer :: ind1,ind2
    integer :: ref_yymmdd,ref_hhmmss
    integer :: it
    integer :: flag
    real(idx) :: sec_start,sec_tmp
    integer :: tmp_yymmdd,tmp_hhmmss
    ind1=1
    ind2=index(time_units,"since")-2
    flag_char=time_units(ind1:ind2)
    ind1=index(time_units,"since")+6
    ind2=ind1+index(time_units(ind1:),"-")-1
    yr_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:),"-")-1
    mn_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:)," ")-1
    dy_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:),":")-1
    hr_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=ind1+index(time_units(ind1:),":")-1
    min_char=time_units(ind1:ind2-1)
    ind1=ind2+1
    ind2=len_trim(time_units)
    sec_char=time_units(ind1:ind2)

    read(yr_char,*) ref_year ; read(mn_char,*) ref_month ; read(dy_char,*) ref_day
    ref_yymmdd=ref_year*10000+ref_month*100+ref_day
    read(hr_char,*) ref_hour ; read(min_char,*) ref_min ; read(sec_char,*) ref_sec
    ref_hhmmss=ref_hour*10000+ref_min*100+ref_sec

    select case(trim(flag_char))
    case("seconds")
       flag=-10000
    case("minitues")
       flag=-100
    case("hours")
       flag=-1
    case("days")
       flag=1
    case("months")
       flag=100
    case("years")
       flag=10000
    end select
    do it = 1,ntime
       call calendar_cal_ymdhms_after(ref_yymmdd,ref_hhmmss,time_in(it),flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,-10000,time(it))
    end do
  end function modify_time
  !==================================================
  ! Read oceanic file
  !==================================================  
  subroutine read_ocean_file(fname,ntime_ocn,time_ocn,temp_in,salt_in,u_in,v_in)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(inout) :: ntime_ocn
    real(idx), allocatable, intent(inout) :: time_ocn(:)
    real(idx),allocatable,intent(inout) :: temp_in(:,:,:,:),salt_in(:,:,:,:),u_in(:,:,:,:),v_in(:,:,:,:)
    integer :: nlon,nlat,nlev
    real(idx), allocatable :: lon_grd(:),lat_grd(:),lev_grd(:)
    integer :: ncid,varid,start(4),count(4)
    !Open file
    call check_r(nf90_open(trim(fname),nf90_nowrite,ncid))
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    ! Obtain level
    call get_dimension(ncid,"lev",nlev,lev_grd)
    ! Obtain time
    call get_dimension(ncid,"time",ntime_ocn,time_ocn)
    ! Obtain T,S,U,V
    start = (/1,1,1,1/)
    count = (/nlon,nlat,nlev,ntime_ocn/)
    allocate(temp_in(nlon,nlat,nlev,ntime_ocn)) ; allocate(salt_in(nlon,nlat,nlev,ntime_ocn))
    allocate(u_in(nlon,nlat,nlev,ntime_ocn)) ; allocate(v_in(nlon,nlat,nlev,ntime_ocn))
    call check_r( nf90_inq_varid(ncid,"temp", varid) )
    call check_r( nf90_get_var(ncid, varid,temp_in, start = start, &
         count = count))
    call check_r( nf90_inq_varid(ncid,"salt", varid) )
    call check_r( nf90_get_var(ncid, varid,salt_in, start = start, &
         count = count))
    call check_r( nf90_inq_varid(ncid,"u", varid) )
    call check_r( nf90_get_var(ncid, varid,u_in, start = start, &
         count = count))
    call check_r( nf90_inq_varid(ncid,"v", varid) )
    call check_r( nf90_get_var(ncid, varid,v_in, start = start, &
         count = count))
    call check_r(nf90_close(ncid))
  end subroutine read_ocean_file
  !==================================================!
  ! Read atmospheric file                            !
  !==================================================!
  subroutine read_atmos_file(fname,varname,ntime,time_atm,data,start_yymmdd,start_hhmmss)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname,varname
    integer,intent(inout) :: ntime
    real(idx), allocatable, intent(inout) :: time_atm(:)
    real(idx),allocatable,intent(inout) :: data(:,:,:)
    integer,intent(in) :: start_yymmdd,start_hhmmss
    integer :: nlon,nlat
    real(idx), allocatable :: lon_grd(:),lat_grd(:),time_grd(:)
    integer :: ncid,varid,start(3),count(3)
    character(len=maxlen) :: time_units
    !Open file
    call check_r( nf90_open(trim(fname),nf90_nowrite,ncid))
    if (allocated(time_grd) .eqv. .true.) then
       deallocate(time_grd)
    end if
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    ! Obtain time
    call get_dimension(ncid,"time",ntime,time_grd)
    call get_var_units(fname,"time",time_units)
    allocate(time_atm(ntime))
    time_atm=modify_time(ntime,time_grd,time_units,start_yymmdd,start_hhmmss)
  ! Obtain atmospheric parameters
    start = (/1,1,1/)
    count = (/nlon,nlat,ntime/)
    if (allocated(data) .eqv. .true.) then
       deallocate(data)
    end if
    allocate(data(nlon,nlat,ntime))
    call check_r(nf90_inq_varid(ncid,trim(varname), varid) )
    call check_r(nf90_get_var(ncid, varid,data, start = start,count = count))
    call check_r(nf90_close(ncid))
  end subroutine read_atmos_file
  subroutine read_atmos_files(fnames,varname,ntime,time_atm,data,start_yymmdd,start_hhmmss)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fnames(:)
    character(len=*),intent(in) :: varname
    integer,intent(inout) :: ntime
    real(idx), allocatable, intent(inout) :: time_atm(:)
    real(idx),allocatable,intent(inout) :: data(:,:,:)
    integer,intent(in) :: start_yymmdd,start_hhmmss
    integer :: nlon,nlat
    real(idx), allocatable :: lon_grd(:),lat_grd(:),time_grd(:),time_tmp(:),data_tmp(:,:,:)
    integer :: ncid,varid,start(3),count(3)
    character(len=maxlen) :: time_units
    integer :: nfile,ifile,ntime_tmp,istr,iend
    nfile=sum(shape(fnames))
    ifile=1
    !Open file
    call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    call check_r(nf90_close(ncid))
    ntime=0
    do ifile=1,nfile
       call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
       call get_dimension(ncid,"time",ntime_tmp,time_grd)
       call get_var_units(fnames(ifile),"time",time_units)
       ntime=ntime+ntime_tmp
       deallocate(time_grd)
       call check_r(nf90_close(ncid))
    end do
    allocate(time_tmp(ntime))
    if (allocated(data) .eqv. .true.) then
       deallocate(data)
    end if
    allocate(data(nlon,nlat,ntime))
    if (allocated(time_atm) .eqv. .true.) then
       deallocate(time_atm)
    end if
    allocate(time_atm(ntime))
    iend=0
    do ifile=1,nfile
       istr=iend+1; 
       call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
       call get_dimension(ncid,"time",ntime_tmp,time_grd)
       start = (/1,1,1/);  count = (/nlon,nlat,ntime_tmp/)
       allocate(data_tmp(1:nlon,1:nlat,1:ntime_tmp))
       iend=iend+ntime_tmp
       time_tmp(istr:iend)=time_grd(1:ntime_tmp)
       call check_r(nf90_inq_varid(ncid,trim(varname), varid) )
       call check_r(nf90_get_var(ncid, varid,data_tmp, start = start,count = count))
       data(1:nlon,1:nlat,istr:iend)=data_tmp(1:nlon,1:nlat,1:ntime_tmp)
       deallocate(time_grd); deallocate(data_tmp)
       call check_r(nf90_close(ncid))
    end do
    time_atm=modify_time(ntime,time_tmp,time_units,start_yymmdd,start_hhmmss)
    deallocate(time_tmp)
  end subroutine read_atmos_files

  !==================================================
  ! Read restart file
  !==================================================  
  subroutine read_init_file(fname,varname,istep,nlon,nlat,nlev,var_in)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname,varname
    integer,intent(in) :: nlon,nlat,nlev,istep
    real(idx),allocatable,intent(inout) :: var_in(:,:,:)
    real(idx), allocatable :: lon_grd(:),lat_grd(:),lev_grd(:)
    integer :: ncid,varid,start(4),count(4)
    !Open file
    call check_r(nf90_open(trim(fname),nf90_nowrite,ncid))
    ! Obtain T,S,U,V
    start = (/1,1,1,istep/)
    count = (/nlon,nlat,nlev,1/)
    if (allocated(var_in) .eqv. .true.) then
       deallocate(var_in)
    end if
    allocate(var_in(nlon,nlat,nlev))
    call check_r( nf90_inq_varid(ncid,varname, varid) )
    call check_r( nf90_get_var(ncid, varid,var_in, start = start, &
         count = count))
    call check_r(nf90_close(ncid))
  end subroutine read_init_file
  !==================================================
  ! Read climatology
  !==================================================  
  subroutine read_clm_file(fname_clm,varname,ntime,time_clm,data_in,start_yymmdd,start_hhmmss)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname_clm,varname
    integer,intent(inout) :: ntime
    real(idx), allocatable, intent(inout) :: time_clm(:)
    real(idx),allocatable,intent(inout) :: data_in(:,:,:,:)
    integer,intent(in) :: start_yymmdd,start_hhmmss
    integer :: nlon,nlat,nlev
    real(idx), allocatable :: lon_grd(:),lat_grd(:),lev_grd(:),time_grd(:)
    integer :: ncid,varid,start(4),count(4)
    character(len=maxlen) :: time_units
    !Open file
    call check_r(nf90_open(trim(fname_clm),nf90_nowrite,ncid))
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    ! Obtain level
    call get_dimension(ncid,"lev",nlev,lev_grd)
    ! Obtain time
    call get_dimension(ncid,"time",ntime,time_grd)
    call get_var_units(fname_clm,"time",time_units)
    allocate(time_clm(ntime))
    time_clm=modify_time(ntime,time_grd,time_units,start_yymmdd,start_hhmmss)
    ! Obtain atmospheric parameters
    ! Obtain T,S,U,V
    start = (/1,1,1,1/)
    count = (/nlon,nlat,nlev,ntime/)
    allocate(data_in(nlon,nlat,nlev,ntime))
    call check_r( nf90_inq_varid(ncid,trim(varname), varid) )
    call check_r( nf90_get_var(ncid, varid,data_in, start = start, &
         count = count))
    call check_r(nf90_close(ncid))
    write(*,*) "*******************************************************"
    write(*,*) " Finish reading "//trim(varname)//" climatological file"
    write(*,*) " File name= "//trim(fname_clm)
    write(*,*) "*******************************************************"
  end subroutine read_clm_file
  subroutine read_3D_files(fnames,varname,ntime,time_atm,data,start_yymmdd,start_hhmmss)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fnames(:)
    character(len=*),intent(in) :: varname
    integer,intent(inout) :: ntime
    real(idx), allocatable, intent(inout) :: time_atm(:)
    real(idx),allocatable,intent(inout) :: data(:,:,:)
    integer,intent(in) :: start_yymmdd,start_hhmmss
    integer :: nlon,nlat
    real(idx), allocatable :: lon_grd(:),lat_grd(:),time_grd(:),time_tmp(:),data_tmp(:,:,:)
    integer :: ncid,varid,start(3),count(3)
    character(len=maxlen) :: time_units
    integer :: nfile,ifile,ntime_tmp,istr,iend
    nfile=sum(shape(fnames))
    ifile=1
    !Open file
    call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    call check_r(nf90_close(ncid))
    ntime=0
    do ifile=1,nfile
       call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
       call get_dimension(ncid,"time",ntime_tmp,time_grd)
       call get_var_units(fnames(ifile),"time",time_units)
       ntime=ntime+ntime_tmp
       deallocate(time_grd)
       call check_r(nf90_close(ncid))
    end do
    allocate(time_tmp(ntime))
    if (allocated(data) .eqv. .true.) then
       deallocate(data)
    end if
    allocate(data(nlon,nlat,ntime))
    if (allocated(time_atm) .eqv. .true.) then
       deallocate(time_atm)
    end if
    allocate(time_atm(ntime))
    iend=0
    do ifile=1,nfile
       istr=iend+1; 
       call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
       call get_dimension(ncid,"time",ntime_tmp,time_grd)
       start = (/1,1,1/);  count = (/nlon,nlat,ntime_tmp/)
       allocate(data_tmp(1:nlon,1:nlat,1:ntime_tmp))
       iend=iend+ntime_tmp
       time_tmp(istr:iend)=time_grd(1:ntime_tmp)
       call check_r(nf90_inq_varid(ncid,trim(varname), varid) )
       call check_r(nf90_get_var(ncid, varid,data_tmp, start = start,count = count))
       data(1:nlon,1:nlat,istr:iend)=data_tmp(1:nlon,1:nlat,1:ntime_tmp)
       deallocate(time_grd); deallocate(data_tmp)
       call check_r(nf90_close(ncid))
       write(*,*) "*******************************************************"
       write(*,*) " Finish reading "//trim(varname)
       write(*,*) " File name= "//trim(fnames(ifile))
       write(*,*) "*******************************************************"
    end do
    time_atm=modify_time(ntime,time_tmp,time_units,start_yymmdd,start_hhmmss)
    deallocate(time_tmp)
  end subroutine read_3D_files

  subroutine read_4D_files(fnames,varname,ntime,time_clm,data,start_yymmdd,start_hhmmss)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fnames(:)
    character(len=*),intent(in) :: varname
    integer,intent(inout) :: ntime
    real(idx), allocatable, intent(inout) :: time_clm(:)
    real(idx),allocatable,intent(inout) :: data(:,:,:,:)
    integer,intent(in) :: start_yymmdd,start_hhmmss
    integer :: nlon,nlat,nlev
    real(idx), allocatable :: lon_grd(:),lat_grd(:),lev_grd(:),time_grd(:),time_tmp(:),data_tmp(:,:,:,:)
    integer :: ncid,varid,start(4),count(4)
    character(len=maxlen) :: time_units
    integer :: nfile,ifile,ntime_tmp,istr,iend
    nfile=sum(shape(fnames))
    ifile=1
    !Open file
    call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    ! Obtain level
    call get_dimension(ncid,"lev",nlev,lev_grd)
    call check_r(nf90_close(ncid))
    ntime=0
    do ifile=1,nfile
       call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
       call get_dimension(ncid,"time",ntime_tmp,time_grd)
       call get_var_units(fnames(ifile),"time",time_units)
       ntime=ntime+ntime_tmp
       deallocate(time_grd)
       call check_r(nf90_close(ncid))
    end do
    allocate(time_tmp(ntime))
    if (allocated(data) .eqv. .true.) then
       deallocate(data)
    end if
    allocate(data(nlon,nlat,nlev,ntime))
    if (allocated(time_clm) .eqv. .true.) then
       deallocate(time_clm)
    end if
    allocate(time_clm(ntime))
    iend=0
    do ifile=1,nfile
       istr=iend+1; 
       call check_r( nf90_open(trim(fnames(ifile)),nf90_nowrite,ncid))
       call get_dimension(ncid,"time",ntime_tmp,time_grd)
       start = (/1,1,1,1/);  count = (/nlon,nlat,nlev,ntime_tmp/)
       allocate(data_tmp(1:nlon,1:nlat,1:nlev,1:ntime_tmp))
       iend=iend+ntime_tmp
       time_tmp(istr:iend)=time_grd(1:ntime_tmp)
       call check_r(nf90_inq_varid(ncid,trim(varname), varid) )
       call check_r(nf90_get_var(ncid, varid,data_tmp, start = start,count = count))
       data(1:nlon,1:nlat,1:nlev,istr:iend)=data_tmp(1:nlon,1:nlat,1:nlev,1:ntime_tmp)
       deallocate(time_grd); deallocate(data_tmp)
       call check_r(nf90_close(ncid))
       write(*,*) "*******************************************************"
       write(*,*) " Finish reading "//trim(varname)
       write(*,*) " File name= "//trim(fnames(ifile))
       write(*,*) "*******************************************************"
    end do
    time_clm=modify_time(ntime,time_tmp,time_units,start_yymmdd,start_hhmmss)
    deallocate(time_tmp)
  end subroutine read_4D_files

  !==================================================
  ! Read climatology
  !==================================================  
  subroutine read_adv_file(fname_adv,varname,ntime,time_adv,data_in,start_yymmdd,start_hhmmss)
    use netcdf
    implicit none
    character(len=*),intent(in) :: fname_adv,varname
    integer,intent(inout) :: ntime
    real(idx), allocatable, intent(inout) :: time_adv(:)
    real(idx),allocatable,intent(inout) :: data_in(:,:,:,:)
    integer,intent(in) :: start_yymmdd,start_hhmmss
    integer :: nlon,nlat,nlev
    real(idx), allocatable :: lon_grd(:),lat_grd(:),lev_grd(:),time_grd(:)
    integer :: ncid,varid,start(4),count(4)
    character(len=maxlen) :: time_units
    !Open file
    call check_r(nf90_open(trim(fname_adv),nf90_nowrite,ncid))
    ! Obtain longitude
    call get_dimension(ncid,"lon",nlon,lon_grd)
    ! Obtain latitude
    call get_dimension(ncid,"lat",nlat,lat_grd)
    ! Obtain level
    call get_dimension(ncid,"lev_rho",nlev,lev_grd)
    ! Obtain time
    call get_dimension(ncid,"time",ntime,time_grd)
    call get_var_units(fname_adv,"time",time_units)
    allocate(time_adv(ntime))
    time_adv=modify_time(ntime,time_grd,time_units,start_yymmdd,start_hhmmss)
    ! Obtain atmospheric parameters
    ! Obtain T,S,U,V
    start = (/1,1,1,1/)
    count = (/nlon,nlat,nlev,ntime/)
    allocate(data_in(nlon,nlat,nlev,ntime))
    call check_r( nf90_inq_varid(ncid,trim(varname), varid) )
    call check_r( nf90_get_var(ncid, varid,data_in, start = start, &
         count = count))
    call check_r(nf90_close(ncid))
    write(*,*) "*******************************************************"
    write(*,*) " Finish reading "//trim(varname)//" correction file"
    write(*,*) " File name= "//trim(fname_adv)
    write(*,*) "*******************************************************"
  end subroutine read_adv_file
  !==================================================================
  ! Linear interpolation(for one value)
  !==================================================================
  subroutine linear_int_one(N,level_in,data_in,newlev,newdata)
    implicit none
    integer :: N
    real(idx),dimension(N),intent(in) :: level_in,data_in
    real(idx),dimension(N) :: level,data
    real(idx) :: newlev,newdata
    integer :: ind_upp,ind_low
    real(idx) :: w1,w2
    integer :: i
    ! Check_R order
    if (level_in(1) .ge. level_in(2)) then
       do i =1,N
          level(i) = level_in(N-i)
          data(i) = data_in(N-i)
       end do
    else
       level = level_in
       data = data_in
    end if
    ! find 2 point that enconpass the value
    if (newlev < minval(level)) then
       ind_low = 1
       ind_upp = 2       
    else if (newlev > maxval(level)) then
       ind_low = N-1
       ind_upp = N
    else
       ind_low=maxloc(level,1,MASK=level<=newlev)
       ind_upp=minloc(level,1,MASK=level>newlev)
    end if
    w1 = newlev - level(ind_low)
    w2 = level(ind_upp) - newlev
    newdata = data(ind_low) * w2 / (w1+w2) +  data(ind_upp) * w1 / (w1+w2)
  end subroutine linear_int_one
  subroutine time_wgt(time_array,time,i1,i2,w1,w2)
    real(idx),intent(in) :: time_array(:),time
    integer,intent(out) :: i1,i2    
    integer :: iasm,ntime
    real(idx) :: w1,w2
    ntime=sum(shape(time_array))
    if (time  .gt. minval(time_array) .and. time .lt. maxval(time_array)) then
       i1 = sum(maxloc(time_array,mask=(time_array<=time)))
       i2 = i1 + 1
       w1 = (time_array(i2)-time) / (time_array(i2)-time_array(i1))
       w2 = (time-time_array(i1)) / (time_array(i2)-time_array(i1))
    else if (time .le. minval(time_array)) then
       i1 = 1
       i2 = 1
       w1 = 1.0_idx
       w2 = 0.0_idx
    else
       i1 = ntime
       i2 = ntime
       w1 = 1.0_idx
       w2 = 0.0_idx
    end if
  end subroutine time_wgt
  subroutine time_wgt2(time_array,time,i1,i2,w1,w2)
    real(idx),intent(in) :: time_array(:),time
    integer,intent(out) :: i1,i2    
    integer :: iasm,ntime
    real(idx) :: w1,w2
    ntime=sum(shape(time_array))
    if (time  .gt. minval(time_array) .and. time .lt. maxval(time_array)) then
       i1 = sum(maxloc(time_array,mask=(time_array<=time)))
       i2 = i1 + 1
       w1 = 0.0_idx
       w2 = 1.0_idx
    else if (time .le. minval(time_array)) then
       i1 = 1
       i2 = 2
       w1 = 1.0_idx
       w2 = 0.0_idx
    else
       i1 = ntime-1
       i2 = ntime
       w1 = 0.0_idx
       w2 = 1.0_idx
    end if
  end subroutine time_wgt2
  function set_data(ind1,ind2,wgt1,wgt2,data_1d) result(data_ret)
    implicit none
    integer,intent(in) :: ind1,ind2
    real(idx),intent(in) :: wgt1,wgt2
    real(idx),intent(in) :: data_1d(:)
    real(idx) :: data_ret
    data_ret= wgt1*data_1d(ind1)+wgt2*data_1d(ind2)
  end function set_data
end module input_files

