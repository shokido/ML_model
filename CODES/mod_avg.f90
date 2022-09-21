module mod_avg
#include 'CPPLISTS.h'
  !================================================
  ! Module program for create/write averaged file
  !================================================
  use param
  use calendar_sub
  use ncdf_write
  implicit none
  ! Time coordinates
  integer :: ntime_avg,iavg
  real(idx),allocatable :: time_avg(:)
  integer,allocatable :: istep_avg(:)
  integer :: out_avg_flag,out_avg_int
  ! rho-file
  character(maxlen) :: fname_out_avg_rho
  real(idx),allocatable :: sw_avg(:,:),lw_avg(:,:),sh_avg(:,:),lh_avg(:,:)
  real(idx),allocatable :: ust_avg(:,:),vst_avg(:,:)
  real(idx),allocatable :: pr_avg(:,:),ev_avg(:,:)
  real(idx),allocatable :: temp_avg(:,:,:),salt_avg(:,:,:),u_avg(:,:,:),v_avg(:,:,:)
  ! q-file
  character(maxlen) :: fname_out_avg_q  
  real(idx),allocatable :: bvf_avg(:,:,:),shear_avg(:,:,:)
  real(idx),allocatable :: qq_avg(:,:,:),l_avg(:,:,:)
  real(idx),allocatable :: km_avg(:,:,:),kt_avg(:,:,:)

  logical :: aout_sw,aout_lw,aout_sh,aout_lh
  logical :: aout_ev,aout_pr,aout_ust,aout_vst
  logical :: aout_temp,aout_salt
  logical :: aout_u,aout_v
  logical :: aout_bvf,aout_shear
  logical :: aout_qq,aout_l
  logical :: aout_km,aout_kt
  
  ! Namelist
  namelist/output_avg_rho_flag/out_avg_flag,out_avg_int
  namelist/output_avg_rho_flag/fname_out_avg_rho
  namelist/output_avg_rho_flag/aout_sw,aout_lw,aout_sh,aout_lh
  namelist/output_avg_rho_flag/aout_ev,aout_pr,aout_ust,aout_vst
  namelist/output_avg_rho_flag/aout_temp,aout_salt,aout_u,aout_v

  namelist/output_avg_q_flag/fname_out_avg_q
  namelist/output_avg_q_flag/aout_bvf,aout_shear
  namelist/output_avg_q_flag/aout_qq,aout_l
  namelist/output_avg_q_flag/aout_km,aout_kt  

  !========================================================
contains
  subroutine create_avg_rho(fname,nlon,nlat,nlev,lon,lat,lev, &
       & start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,out_int,&
       & out_sw,out_lw,out_sh,out_lh,&
       & out_ev,out_pr,out_ust,out_vst,&
       & out_temp,out_salt,out_u,out_v,&
       & missing_value,istep,ntime,time)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev
    real(idx),intent(in) :: lon(nlon),lat(nlat),lev(nlev)
    integer,intent(in) :: start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss
    integer,intent(in) :: out_flag,out_int
    logical,intent(in) :: out_temp,out_salt,out_u,out_v
    logical,intent(in) :: out_sw,out_lw,out_sh,out_lh
    logical,intent(in) :: out_ev,out_pr,out_ust,out_vst
    real(idx),intent(in) :: missing_value
    integer,allocatable,intent(inout) :: istep(:)
    integer,intent(inout) :: ntime
    real(idx),allocatable,intent(inout) :: time(:)
    character(len=maxlen) :: ref_time
    real(idx) :: tmp
    integer :: it
    integer :: tmp_yymmdd,tmp_hhmmss
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,out_flag)
    ref_time=calendar_create_time_att(start_yymmdd,start_hhmmss,out_flag)
    call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,end_yymmdd,end_hhmmss,out_flag,tmp)
    ntime=int(tmp / out_int)
    ntime=max(ntime,1)
    allocate(time(ntime)); allocate(istep(ntime))
    do it=1,ntime
       time(it) = real(it)* out_int
       call calendar_cal_ymdhms_after(start_yymmdd,start_hhmmss,time(it),out_flag,tmp_yymmdd,tmp_hhmmss)
       call calendar_cal_length_ymdhms(start_yymmdd,start_hhmmss,tmp_yymmdd,tmp_hhmmss,-10000,tmp)
       istep(it)=int(tmp/dt)
    end do
    write(*,*) "Create avg_rho=",trim(fname)
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","lev","time",&
         & "degrees_east","degrees_north","m",trim(ref_time), &
         & lon,lat,lev,time)
    if (out_sw .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"sw"/),(/"W/m^2"/),missing_value)
    end if
    if (out_lw .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"lw"/),(/"W/m^2"/),missing_value)
    end if
   if (out_sh .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"sh"/),(/"W/m^2"/),missing_value)
    end if
    if (out_lh .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"lh"/),(/"W/m^2"/),missing_value)
    end if
    if (out_ev .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"ev"/),(/"m/s"/),missing_value)
    end if
    if (out_pr .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"pr"/),(/"m/s"/),missing_value)
    end if
    if (out_ust .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"ust"/),(/"N/m^2"/),missing_value)
    end if
    if (out_vst .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","time",&
            & 1,(/"vst"/),(/"N/m^2"/),missing_value)
    end if
    if (out_temp .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"temp"/),(/"degrees celcius"/),missing_value)
    end if
    if (out_salt .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"salt"/),(/"psu"/),missing_value)
    end if
    if (out_u .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"u"/),(/"m/s"/),missing_value)
    end if
    if (out_v .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"v"/),(/"m/s"/),missing_value)
    end if
  end subroutine create_avg_rho
  ! Output q-lev file
  subroutine create_avg_q(fname,nlon,nlat,nlev,ntime,lon,lat,z_q,time, &
       & start_yymmdd,start_hhmmss,flag_yymmdd,out_bvf,out_shear,out_qq,out_l,out_km,out_kt,missing_value)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,ntime
    real(idx),intent(in) :: lon(nlon),lat(nlat),z_q(nlev),time(ntime)
    integer,intent(in) :: start_yymmdd,start_hhmmss,flag_yymmdd
    logical,intent(in) :: out_bvf,out_shear,out_qq,out_l
    logical,intent(in) :: out_km,out_kt
    real(idx) :: missing_value
    character :: time_unit*100
    write(*,*) "Create avg_rho=",trim(fname)
    time_unit=calendar_create_time_att(start_yymmdd,start_hhmmss,flag_yymmdd)
    call writenet_pre(fname,nlon,nlat,nlev,ntime,"lon","lat","lev","time",&
         & "degrees_east","degrees_north","m",trim(time_unit), &
         & lon,lat,z_q,time)
    if (out_bvf .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"bvf"/),(/"1/s^2"/),missing_value)
    end if
    if (out_shear .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"shear"/),(/"1/s^2"/),missing_value)
    end if
#if defined(NNF) | defined(KC) | defined(bd)
    if (out_qq .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"qq"/),(/"m^2/s^2"/),missing_value)
    end if
    if (out_l .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"l"/),(/"m"/),missing_value)
    end if
#endif
    if (out_km .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"km"/),(/"m^2/s"/),missing_value)
    end if
    if (out_kt .eqv. .TRUE.) then
       call writenet_dv(trim(fname),"lon","lat","lev","time",&
            & 1,(/"kt"/),(/"m^2/s"/),missing_value)
    end if
  end subroutine create_avg_q

  subroutine output_avg_rho(fname,nlon,nlat,nlev,itime,&
       & out_sw,out_lw,out_sh,out_lh,out_ev,out_pr,&
       & out_ust,out_vst,out_temp,out_salt,out_u,out_v,&
       & sw,lw,sh,lh,ev,pr,ust,vst,&
       & temp,salt,u,v,avg_count)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_sw,out_lw,out_sh,out_lh
    logical,intent(in) :: out_ev,out_pr,out_ust,out_vst
    logical,intent(in) :: out_temp,out_salt,out_u,out_v
    real(idx),intent(inout) :: sw(nlon,nlat),lw(nlon,nlat),sh(nlon,nlat),lh(nlon,nlat)
    real(idx),intent(inout) :: ev(nlon,nlat),pr(nlon,nlat),ust(nlon,nlat),vst(nlon,nlat)
    real(idx),intent(inout) :: temp(nlon,nlat,nlev),salt(nlon,nlat,nlev),u(nlon,nlat,nlev),v(nlon,nlat,nlev)
    integer,intent(in) :: avg_count
    real(idx) :: tmp_3D(1:nlon,1:nlat,1)
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    ! Divide
    sw = sw / avg_count ; lw = lw / avg_count
    sh = sh / avg_count ; lh = lh / avg_count
    ev = ev / avg_count ; pr = pr / avg_count
    ust = ust / avg_count ; vst = vst / avg_count
    temp =  temp / avg_count ; salt =  salt / avg_count
    u = u / avg_count ; v =  v / avg_count 
    tmp_3D = 0.0_idx
    tmp_4D = 0.0_idx
    if (out_sw .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=sw(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"sw",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_lw .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=lw(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"lw",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_sh .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=sh(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"sh",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_lh .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=lh(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"lh",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_ev .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=ev(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"ev",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_pr .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=pr(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"pr",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_ust .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=ust(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"ust",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if
    if (out_vst .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=vst(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"vst",1,nlon,1,nlat,itime,itime,&
            & tmp_3D)
    end if

    if (out_temp .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"temp",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    if (out_salt .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=salt(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"salt",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    if (out_u .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=u(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"u",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    if (out_v .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=v(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"v",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    sw = 0.0_idx ; lw = 0.0_idx ; sh = 0.0_idx ; lh = 0.0_idx
    ev = 0.0_idx ; pr = 0.0_idx ; ust = 0.0_idx ; vst = 0.0_idx
    temp =  0.0_idx ; salt =  0.0_idx
    u = 0.0_idx ; v = 0.0_idx
  end subroutine output_avg_rho
  ! Output
  subroutine output_avg_q(fname,nlon,nlat,nlev,itime,&
       & out_bvf,out_shear,out_qq,out_l,out_km,out_kt,bvf,shear,&
#if defined(NNF) | defined(KC) | defined(bd)
       & qq,l,&
#endif
       & km,kt,avg_count)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_bvf,out_shear,out_qq,out_l
    logical,intent(in) :: out_km,out_kt
    real(idx),intent(inout) :: bvf(nlon,nlat,nlev),shear(nlon,nlat,nlev)
#if defined(NNF) | defined(KC) | defined(bd)
    real(idx),intent(inout) :: qq(nlon,nlat,nlev),l(nlon,nlat,nlev)
#endif
    real(idx),intent(inout) :: km(nlon,nlat,nlev),kt(nlon,nlat,nlev)
    integer,intent(in) :: avg_count
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    ! Divide
    bvf =  bvf / avg_count ; shear =  shear / avg_count
#if defined(NNF) | defined(KC) | defined(bd)
    qq =  qq / avg_count ; l =  l / avg_count
#endif
    km =  km / avg_count ; kt =  kt / avg_count

    tmp_4D=0.0_idx
    if (out_bvf .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=bvf(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"bvf",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    if (out_shear .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=shear(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"shear",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
#if defined(NNF) | defined(KC) | defined(bd)
    if (out_qq .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    if (out_l .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=l(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"l",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
#endif
    if (out_km .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=km(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"km",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    if (out_kt .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=kt(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"kt",1,nlon,1,nlat,1,nlev,itime,itime,&
            & tmp_4D)
    end if
    bvf = 0.0_idx ; shear = 0.0_idx
#if defined(NNF) | defined(KC) | defined(bd)
    qq = 0.0_idx ; l = 0.0_idx
#endif
    km =  0.0_idx ; kt = 0.0_idx
  end subroutine output_avg_q
end module mod_avg
