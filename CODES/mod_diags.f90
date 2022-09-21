module mod_diags
  !==================================================!
  ! Module program for create/write diagnositic file !
  !==================================================!
  use param
  use calendar_sub
  use ncdf_write
  implicit none
  ! Temp
  integer :: ntime_tdia,itdia
  real(idx),allocatable :: time_tdia(:)
  integer,allocatable :: istep_tdia(:)
  integer :: out_tdia_flag,out_tdia_int
  character(maxlen) :: fname_out_tdiag
  real(idx),allocatable :: temp_rate(:,:,:),temp_vdiff(:,:,:)
  real(idx),allocatable :: temp_penet(:,:,:),temp_relax(:,:,:)
  logical :: dout_temp_rate,dout_temp_vdiff
  logical :: dout_temp_penet,dout_temp_relax

  ! Salinity
  integer :: ntime_sdia,isdia
  real(idx),allocatable :: time_sdia(:)
  integer,allocatable :: istep_sdia(:)
  integer :: out_sdia_flag,out_sdia_int
  character(maxlen) :: fname_out_sdiag
  real(idx),allocatable :: salt_rate(:,:,:),salt_vdiff(:,:,:),salt_relax(:,:,:)
  logical :: dout_salt_rate,dout_salt_vdiff
  logical :: dout_salt_relax

  ! Velocity
  integer :: ntime_uvdia,iuvdia
  real(idx),allocatable :: time_uvdia(:)
  integer,allocatable :: istep_uvdia(:)
  integer :: out_uvdia_flag,out_uvdia_int
  character(maxlen) :: fname_out_uvdiag
  real(idx),allocatable :: u_rate(:,:,:),u_cor(:,:,:),u_vdiff(:,:,:),u_relax(:,:,:)
  real(idx),allocatable :: v_rate(:,:,:),v_cor(:,:,:),v_vdiff(:,:,:),v_relax(:,:,:)
  logical :: dout_uv_rate,dout_uv_cor,dout_uv_vdiff,dout_uv_relax

  ! TKE
  integer :: ntime_qqdia,iqqdia
  real(idx),allocatable :: time_qqdia(:)
  integer,allocatable :: istep_qqdia(:)
  integer :: out_qqdia_flag,out_qqdia_int
  character(maxlen) :: fname_out_qqdiag
  real(idx),allocatable :: qq_rate(:,:,:),qq_vdiff(:,:,:)
  real(idx),allocatable :: qq_sp(:,:,:),qq_bp(:,:,:),qq_disp(:,:,:)
  logical :: dout_qq_rate,dout_qq_vdiff
  logical :: dout_qq_sp,dout_qq_bp
  logical :: dout_qq_disp  

  ! Namelist
  namelist/output_temp_diag_flag/out_tdia_flag,out_tdia_int
  namelist/output_temp_diag_flag/fname_out_tdiag
  namelist/output_temp_diag_flag/dout_temp_rate,dout_temp_vdiff
  namelist/output_temp_diag_flag/dout_temp_penet,dout_temp_relax

  namelist/output_salt_diag_flag/out_sdia_flag,out_sdia_int
  namelist/output_salt_diag_flag/fname_out_sdiag
  namelist/output_salt_diag_flag/dout_salt_rate,dout_salt_vdiff,dout_salt_relax

  namelist/output_uv_diag_flag/out_uvdia_flag,out_uvdia_int
  namelist/output_uv_diag_flag/fname_out_uvdiag
  namelist/output_uv_diag_flag/dout_uv_rate,dout_uv_cor,dout_uv_vdiff,dout_uv_relax
  namelist/output_qq_diag_flag/out_qqdia_flag,out_qqdia_int
  namelist/output_qq_diag_flag/fname_out_qqdiag
  namelist/output_qq_diag_flag/dout_qq_rate,dout_qq_vdiff
  namelist/output_qq_diag_flag/dout_qq_sp,dout_qq_bp
  namelist/output_qq_diag_flag/dout_qq_disp

contains
  ! Allocate arrays
  subroutine allocate_tdiag(nlon,nlat,nlev,temp_rate,temp_vdiff,temp_penet,temp_relax)
    implicit none
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),allocatable,intent(inout) :: temp_rate(:,:,:),temp_vdiff(:,:,:)
    real(idx),allocatable,intent(inout) :: temp_penet(:,:,:),temp_relax(:,:,:)
    allocate(temp_rate(nlon,nlat,nlev))
    allocate(temp_vdiff(nlon,nlat,nlev))
    allocate(temp_penet(nlon,nlat,nlev))
    allocate(temp_relax(nlon,nlat,nlev))
    temp_rate = 0.0_idx
    temp_vdiff = 0.0_idx
    temp_penet = 0.0_idx
    temp_relax = 0.0_idx
  end subroutine allocate_tdiag
  subroutine deallocate_tdiag(temp_rate,temp_vdiff,temp_penet,temp_relax)
    implicit none
    real(idx),allocatable,intent(inout) :: temp_rate(:,:,:),temp_vdiff(:,:,:)
    real(idx),allocatable,intent(inout) :: temp_penet(:,:,:),temp_relax(:,:,:)
    deallocate(temp_rate); deallocate(temp_vdiff)
    deallocate(temp_penet); deallocate(temp_relax)
  end subroutine deallocate_tdiag
  subroutine allocate_sdiag(nlon,nlat,nlev,salt_rate,salt_vdiff,salt_relax)
    implicit none
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),allocatable,intent(inout) :: salt_rate(:,:,:),salt_vdiff(:,:,:),salt_relax(:,:,:)
    allocate(salt_rate(nlon,nlat,nlev))
    allocate(salt_vdiff(nlon,nlat,nlev))
    allocate(salt_relax(nlon,nlat,nlev))
    salt_rate = 0.0_idx
    salt_vdiff = 0.0_idx
    salt_relax = 0.0_idx
  end subroutine allocate_sdiag
  subroutine deallocate_sdiag(salt_rate,salt_vdiff,salt_relax)
    implicit none
    real(idx),allocatable,intent(inout) :: salt_rate(:,:,:),salt_vdiff(:,:,:),salt_relax(:,:,:)
    deallocate(salt_rate)
    deallocate(salt_vdiff)
    deallocate(salt_relax)
  end subroutine deallocate_sdiag
  subroutine allocate_uvdiag(nlon,nlat,nlev,u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
    implicit none
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),allocatable,intent(inout) :: u_rate(:,:,:),u_cor(:,:,:),u_vdiff(:,:,:),u_relax(:,:,:)
    real(idx),allocatable,intent(inout) :: v_rate(:,:,:),v_cor(:,:,:),v_vdiff(:,:,:),v_relax(:,:,:)
    allocate(u_rate(nlon,nlat,nlev)); allocate(u_cor(nlon,nlat,nlev))
    allocate(u_vdiff(nlon,nlat,nlev)); allocate(u_relax(nlon,nlat,nlev))
    allocate(v_rate(nlon,nlat,nlev)); allocate(v_cor(nlon,nlat,nlev))
    allocate(v_vdiff(nlon,nlat,nlev)); allocate(v_relax(nlon,nlat,nlev))
    u_rate = 0.0_idx ; u_cor = 0.0_idx ; u_vdiff = 0.0_idx ; u_relax = 0.0_idx
    v_rate = 0.0_idx ; v_cor = 0.0_idx ; v_vdiff = 0.0_idx ; v_relax = 0.0_idx
  end subroutine allocate_uvdiag
  subroutine deallocate_uvdiag(u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
    implicit none
    real(idx),allocatable,intent(inout) :: u_rate(:,:,:),u_cor(:,:,:),u_vdiff(:,:,:),u_relax(:,:,:)
    real(idx),allocatable,intent(inout) :: v_rate(:,:,:),v_cor(:,:,:),v_vdiff(:,:,:),v_relax(:,:,:)
    deallocate(u_rate); deallocate(u_cor)
    deallocate(u_vdiff); deallocate(u_relax)
    deallocate(v_rate); deallocate(v_cor)
    deallocate(v_vdiff); deallocate(v_relax)
  end subroutine deallocate_uvdiag
  subroutine allocate_qqdiag(nlon,nlat,nlev,qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp)
    implicit none
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),allocatable,intent(inout) :: qq_rate(:,:,:),qq_vdiff(:,:,:),qq_sp(:,:,:),qq_bp(:,:,:),qq_disp(:,:,:)
    allocate(qq_rate(nlon,nlat,nlev))
    allocate(qq_vdiff(nlon,nlat,nlev))
    allocate(qq_sp(nlon,nlat,nlev))
    allocate(qq_bp(nlon,nlat,nlev))
    allocate(qq_disp(nlon,nlat,nlev))
    qq_rate = 0.0_idx ;  qq_vdiff = 0.0_idx ; qq_sp = 0.0_idx ; qq_bp=0.0_idx ; qq_disp = 0.0_idx
  end subroutine allocate_qqdiag
  subroutine deallocate_qqdiag(qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp)
    implicit none
    real(idx),allocatable,intent(inout) :: qq_rate(:,:,:),qq_vdiff(:,:,:),qq_sp(:,:,:),qq_bp(:,:,:),qq_disp(:,:,:)
    deallocate(qq_rate)
    deallocate(qq_vdiff)
    deallocate(qq_sp)
    deallocate(qq_bp)
    deallocate(qq_disp)
  end subroutine deallocate_qqdiag
  subroutine create_tdiag(fname,nlon,nlat,nlev,lon,lat,lev, &
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss, &
       & out_flag,out_int,out_trate,out_tvdiff,out_tpenet,out_trelax,&
       & missing_value,istep,ntime,time)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev)
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
    integer,intent(in) :: out_flag,out_int
    logical,intent(in) :: out_trate
    logical,intent(in) :: out_tvdiff
    logical,intent(in) :: out_tpenet
    logical,intent(in) :: out_trelax
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: ntime
    real(idx),allocatable,intent(inout) :: time(:)
    character(len=maxlen) :: ref_time
    real(idx) :: tmp
    integer :: it
    integer :: tmp_yymmdd,tmp_hhmmss
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,out_flag)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    ntime= int(tmp / out_int)
    ntime = max(ntime,1)
    allocate(time(ntime)); allocate(istep(ntime))
    do it=1,ntime
       time(it) = (real(it))* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,-10000,tmp)
       istep(it)=int(tmp/dt)
    end do
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","lev","time",&
         & "degrees_east","degrees_north","m",trim(ref_time), &
         & lon,lat,lev,time)
    if (out_trate .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"temp_rate"/),(/"K/s"/),missing_value)
    end if
    if (out_tvdiff .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"temp_vdiff"/),(/"K/s"/),missing_value)
    end if
    if (out_tpenet .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"temp_penet"/),(/"K/s"/),missing_value)
    end if
    if (out_trelax .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"temp_relax"/),(/"K/s"/),missing_value)
    end if
  end subroutine create_tdiag
  ! Salinity diagnostics
  subroutine create_sdiag(fname,nlon,nlat,nlev,lon,lat,lev, &
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss, &
       & out_flag,out_int,out_srate,out_svdiff,out_srelax,&
       & missing_value,istep,ntime,time)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev)
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
    integer,intent(in) :: out_flag,out_int
    logical,intent(in) :: out_srate
    logical,intent(in) :: out_svdiff
    logical,intent(in) :: out_srelax
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: ntime
    real(idx),allocatable,intent(inout) :: time(:)
    character(len=maxlen) :: ref_time
    real(idx) :: tmp
    integer :: it
    integer :: tmp_yymmdd,tmp_hhmmss
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,out_flag)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    ntime= int(tmp / out_int)
    ntime = max(ntime,1)
    allocate(time(ntime)); allocate(istep(ntime))
    do it=1,ntime
       time(it) = real(it)* out_int       
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,-10000,tmp)
       istep(it)=int(tmp/dt)
    end do
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","lev","time",&
         & "degrees_east","degrees_north","m",trim(ref_time), &
         & lon,lat,lev,time)
    if (out_srate .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"salt_rate"/),(/"psu/s"/),missing_value)
    end if
    if (out_svdiff .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"salt_vdiff"/),(/"psu/s"/),missing_value)
    end if
    if (out_srelax .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"salt_relax"/),(/"psu/s"/),missing_value)
    end if
  end subroutine create_sdiag
  ! Velocity diagnostics
  subroutine create_uvdiag(fname,nlon,nlat,nlev,lon,lat,lev, &
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss, &
       & out_flag,out_int,out_uvrate,out_uvcor,out_uvvdiff,out_uvrelax,&
       & missing_value,istep,ntime,time)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev)
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
    integer,intent(in) :: out_flag,out_int
    logical,intent(in) :: out_uvrate,out_uvcor,out_uvvdiff,out_uvrelax
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: ntime
    real(idx),allocatable,intent(inout) :: time(:)
    character(len=maxlen) :: ref_time
    real(idx) :: tmp
    integer :: it
    integer :: tmp_yymmdd,tmp_hhmmss
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,out_flag)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    ntime= int(tmp / out_int)
    ntime = max(ntime,1)
    allocate(time(ntime)); allocate(istep(ntime))
    do it=1,ntime
       time(it) = real(it)* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,-10000,tmp)
       istep(it)=int(tmp/dt)
    end do
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","lev","time",&
         & "degrees_east","degrees_north","m",trim(ref_time), &
         & lon,lat,lev,time)
    if (out_uvrate .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 2,(/"u_rate","v_rate"/),(/"m/s^2","m/s^2"/),missing_value)
    end if
    if (out_uvcor .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 2,(/"u_cor","v_cor"/),(/"m/s^2","m/s^2"/),missing_value)
    end if
    if (out_uvvdiff .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 2,(/"u_vdiff","v_vdiff"/),(/"m/s^2","m/s^2"/),missing_value)
    end if
    if (out_uvrelax .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 2,(/"u_relax","v_relax"/),(/"m/s^2","m/s^2"/),missing_value)
    end if
  end subroutine create_uvdiag
  subroutine create_qqdiag(fname,nlon,nlat,nlev,lon,lat,lev, &
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss, &
       & out_flag,out_int,out_qqrate,out_qqvdiff,out_qqsp,out_qqbp,out_qqdisp,&
       & missing_value,istep,ntime,time)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev)
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
    integer,intent(in) :: out_flag,out_int
    logical,intent(in) :: out_qqrate,out_qqvdiff,out_qqsp,out_qqbp,out_qqdisp
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: ntime
    real(idx),allocatable,intent(inout) :: time(:)
    character(len=maxlen) :: ref_time
    real(idx) :: tmp
    integer :: it
    integer :: tmp_yymmdd,tmp_hhmmss
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,out_flag)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    ntime= int(tmp / out_int)
    ntime = max(ntime,1)
    allocate(time(ntime)); allocate(istep(ntime))
    do it=1,ntime
       time(it) = real(it)* out_int       
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,-10000,tmp)
       istep(it)=int(tmp/dt)
    end do
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","lev","time",&
         & "degrees_east","degrees_north","m",trim(ref_time), &
         & lon,lat,lev,time)
    if (out_qqrate .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"qq_rate"/),(/"m^2/s^3"/),missing_value)
    end if
    if (out_qqvdiff .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"qq_vdiff"/),(/"m^2/s^3"/),missing_value)
    end if
    if (out_qqsp .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"qq_sp"/),(/"m^2/s^3"/),missing_value)
    end if
    if (out_qqbp .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"qq_bp"/),(/"m^2/s^3"/),missing_value)
    end if
    if (out_qqdisp .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"qq_disp"/),(/"m^2/s^3"/),missing_value)
    end if

  end subroutine create_qqdiag

  !===========================================================
  ! Output temp-diag
  !===========================================================
  subroutine output_diag_temp(fname,nlon,nlat,nlev,itime,&
       & out_trate,out_tvdiff,out_tpenet,out_trelax,&
       & temp_rate,temp_vdiff,temp_penet,temp_relax,diag_count)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_trate
    logical,intent(in) :: out_tvdiff
    logical,intent(in) :: out_tpenet
    logical,intent(in) :: out_trelax
    real(idx),intent(inout) :: temp_rate(nlon,nlat,nlev),temp_vdiff(nlon,nlat,nlev)
    real(idx),intent(inout) :: temp_penet(nlon,nlat,nlev),temp_relax(nlon,nlat,nlev)
    integer,intent(in) :: diag_count
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    ! Divide
    temp_rate=temp_rate / diag_count ; temp_vdiff=temp_vdiff / diag_count
    temp_penet=temp_penet / diag_count ;  temp_relax=temp_relax / diag_count
    
    tmp_4D=0.0_idx
    if (out_trate .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp_rate(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"temp_rate",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_tvdiff .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp_vdiff(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"temp_vdiff",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_tpenet .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp_penet(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"temp_penet",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_trelax .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp_relax(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"temp_relax",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    temp_rate(:,:,:)=0.0_idx ; temp_vdiff(:,:,:)=0.0_idx
    temp_penet(:,:,:)=0.0_idx ; temp_relax(:,:,:)=0.0_idx
  end subroutine output_diag_temp
  !===========================================================
  ! Output salt-diag
  !===========================================================
  subroutine output_diag_salt(fname,nlon,nlat,nlev,itime,&
       & out_srate,out_svdiff,out_srelax,&
       & salt_rate,salt_vdiff,salt_relax,diag_count)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_srate
    logical,intent(in) :: out_svdiff
    logical,intent(in) :: out_srelax
!    real(idx),intent(in) :: 
    real(idx),intent(inout) :: salt_rate(nlon,nlat,nlev),salt_vdiff(nlon,nlat,nlev)
    real(idx),intent(inout) :: salt_relax(nlon,nlat,nlev)
    integer,intent(in) :: diag_count
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    ! Divide
    salt_rate=salt_rate / diag_count ; salt_vdiff=salt_vdiff / diag_count
    salt_relax=salt_relax / diag_count

    tmp_4D=0.0_idx
    if (out_srate .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=salt_rate(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"salt_rate",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_svdiff .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=salt_vdiff(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"salt_vdiff",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_srelax .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=salt_relax(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"salt_relax",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    salt_rate(:,:,:)=0.0_idx ; salt_vdiff(:,:,:)=0.0_idx
    salt_relax(:,:,:)=0.0_idx
  end subroutine output_diag_salt
  !===========================================================
  ! Output uv-diag
  !===========================================================
  subroutine output_diag_uv(fname,nlon,nlat,nlev,itime,&
       & out_uvrate,out_uvcor,out_uvvdiff,out_uvrelax,&
       & u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax,diag_count)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_uvrate,out_uvcor,out_uvvdiff,out_uvrelax
    real(idx),intent(inout) :: u_rate(nlon,nlat,nlev),u_cor(nlon,nlat,nlev)
    real(idx),intent(inout) :: u_vdiff(nlon,nlat,nlev),u_relax(nlon,nlat,nlev)
    real(idx),intent(inout) :: v_rate(nlon,nlat,nlev),v_cor(nlon,nlat,nlev)
    real(idx),intent(inout) :: v_vdiff(nlon,nlat,nlev),v_relax(nlon,nlat,nlev)
    integer,intent(in) :: diag_count
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    ! Divide
    u_rate=u_rate / diag_count ;  u_cor=u_cor / diag_count
    u_vdiff=u_vdiff / diag_count; u_relax=u_relax / diag_count
    v_rate=v_rate / diag_count ;  v_cor=v_cor / diag_count
    v_vdiff=v_vdiff / diag_count; v_relax=v_relax / diag_count

    tmp_4D=0.0_idx
    if (out_uvrate .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=u_rate(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"u_rate",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=v_rate(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"v_rate",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_uvcor .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=u_cor(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"u_cor",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=v_cor(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"v_cor",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_uvvdiff .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=u_vdiff(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"u_vdiff",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=v_vdiff(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"v_vdiff",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_uvrelax .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=u_relax(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"u_relax",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=v_relax(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"v_relax",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    u_rate(:,:,:)=0.0_idx ; u_cor(:,:,:)=0.0_idx
    u_vdiff(:,:,:)=0.0_idx; u_relax(:,:,:)=0.0
    v_rate(:,:,:)=0.0_idx ; v_cor(:,:,:)=0.0_idx
    v_vdiff(:,:,:)=0.0_idx; v_relax(:,:,:)=0.0
  end subroutine output_diag_uv
  !===========================================================
  ! Output salt-diag
  !===========================================================
  subroutine output_diag_qq(fname,nlon,nlat,nlev,itime,&
       & out_qqrate,out_qqvdiff,out_qqsp,out_qqbp,out_qqdisp,&
       & qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp,diag_count)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_qqrate,out_qqvdiff,out_qqsp,out_qqbp,out_qqdisp
!    real(idx),intent(in) :: 
    real(idx),intent(inout) :: qq_rate(nlon,nlat,nlev),qq_vdiff(nlon,nlat,nlev)
    real(idx),intent(inout) :: qq_sp(nlon,nlat,nlev),qq_bp(nlon,nlat,nlev),qq_disp(nlon,nlat,nlev)
    integer,intent(in) :: diag_count
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    ! Divide
    qq_rate=qq_rate / diag_count ; qq_vdiff=qq_vdiff / diag_count
    qq_sp=qq_sp / diag_count; qq_bp=qq_bp / diag_count; qq_disp = qq_disp / diag_count

    tmp_4D=0.0_idx
    if (out_qqrate .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq_rate(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq_rate",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_qqvdiff .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq_vdiff(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq_vdiff",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_qqsp .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq_sp(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq_sp",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_qqbp .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq_bp(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq_bp",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_qqdisp .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq_disp(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq_disp",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    qq_rate(:,:,:)=0.0_idx ; qq_vdiff(:,:,:)=0.0_idx
    qq_sp(:,:,:)=0.0_idx ; qq_bp(:,:,:)=0.0_idx
    qq_disp(:,:,:)=0.0_idx 
  end subroutine output_diag_qq
end module mod_diags
