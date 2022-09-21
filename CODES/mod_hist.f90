module mod_hist
#include "CPPLISTS.h"
  !================================================
  ! Module program for create/write history file
  ! (Snapshot file)
  !================================================
  use param
  use calendar_sub
  use ncdf_write
  implicit none
  ! Time coordinates
  integer :: ntime_hist,ihist
  real(idx),allocatable :: time_hist(:)
  integer,allocatable :: istep_hist(:)
  integer :: out_hist_flag,out_hist_int
  ! rho-file
  character(maxlen) :: fname_out_hist_rho
  real(idx),allocatable :: sw(:,:),lw(:,:),sh(:,:),lh(:,:)
  real(idx),allocatable :: ust(:,:),vst(:,:)
  real(idx),allocatable :: pr(:,:),ev(:,:)
  real(idx),allocatable :: temp(:,:,:),salt(:,:,:),u(:,:,:),v(:,:,:)
  real(idx),allocatable :: temp_past(:,:,:),salt_past(:,:,:),u_past(:,:,:),v_past(:,:,:)
  ! q-file
  character(maxlen) :: fname_out_hist_q
  real(idx),allocatable :: bvf(:,:,:),shear(:,:,:)
  real(idx),allocatable :: qq(:,:,:),l(:,:,:)
  real(idx),allocatable :: km(:,:,:),kt(:,:,:)
  real(idx),allocatable :: qq_past(:,:,:)
  
  ! Output at z_rho
  logical :: hout_sw,hout_lw,hout_sh,hout_lh
  logical :: hout_ev,hout_pr,hout_ust,hout_vst
  logical :: hout_temp,hout_salt
  logical :: hout_u,hout_v

  ! Output at z_q
  logical :: hout_bvf,hout_shear
  logical :: hout_qq,hout_l
  logical :: hout_km,hout_kt

  ! Namelist
  namelist/output_hist_rho_flag/out_hist_flag,out_hist_int 
  namelist/output_hist_rho_flag/fname_out_hist_rho
  namelist/output_hist_rho_flag/hout_sw,hout_lw,hout_sh,hout_lh
  namelist/output_hist_rho_flag/hout_ev,hout_pr,hout_ust,hout_vst
  namelist/output_hist_rho_flag/hout_temp,hout_salt
  namelist/output_hist_rho_flag/hout_u,hout_v
  namelist/output_hist_q_flag/fname_out_hist_q
  namelist/output_hist_q_flag/hout_bvf,hout_shear
  namelist/output_hist_q_flag/hout_qq,hout_l
  namelist/output_hist_q_flag/hout_km,hout_kt
  
  !========================================================
contains
  ! Allocate arrays
  subroutine allocate_ocn_3d_arrays(nlon,nlat,nz,temp,salt,u,v,qq,l)
    implicit none
    integer,intent(in) :: nlon,nlat,nz
    real(idx),allocatable :: temp(:,:,:),salt(:,:,:),u(:,:,:),v(:,:,:)
    real(idx),allocatable :: qq(:,:,:),l(:,:,:)
    allocate(temp(nlon,nlat,nz)) ; allocate(salt(nlon,nlat,nz))
    allocate(u(nlon,nlat,nz)) ;  allocate(v(nlon,nlat,nz))
    allocate(qq(nlon,nlat,nz)) ;  allocate(l(nlon,nlat,nz))
    temp = 0.0_idx ; salt =  0.0_idx ; u= 0.0_idx ; v= 0.0_idx
    qq= 0.0_idx ; l= 0.0_idx
  end subroutine allocate_ocn_3d_arrays
  subroutine deallocate_ocn_3d_arrays(temp,salt,u,v,qq,l)
    implicit none
    real(idx),allocatable :: temp(:,:,:),salt(:,:,:),u(:,:,:),v(:,:,:)
    real(idx),allocatable :: qq(:,:,:),l(:,:,:)
    deallocate(temp) ; deallocate(salt)
    deallocate(u) ;  deallocate(v)
    deallocate(qq) ;  deallocate(l)
  end subroutine deallocate_ocn_3d_arrays
  subroutine allocate_atm_2d_arrays(nlon,nlat,sw,lw,sh,lh,ev,pr,ust,vst)
    implicit none
    integer,intent(in) :: nlon,nlat
    real(idx),allocatable :: sw(:,:),lw(:,:),sh(:,:),lh(:,:)
    real(idx),allocatable :: ev(:,:),pr(:,:)
    real(idx),allocatable :: ust(:,:),vst(:,:)
    allocate(sw(nlon,nlat)); allocate(lw(nlon,nlat))
    allocate(sh(nlon,nlat)); allocate(lh(nlon,nlat))
    allocate(ev(nlon,nlat)); allocate(pr(nlon,nlat))
    allocate(ust(nlon,nlat)); allocate(vst(nlon,nlat))
    sw=0.0_idx ; lw=0.0_idx ; sh=0.0_idx ; lh=0.0_idx
    ev=0.0_idx ; pr=0.0_idx
    ust=0.0_idx ; vst=0.0_idx
  end subroutine allocate_atm_2d_arrays
  subroutine deallocate_atm_2d_arrays(sw,lw,sh,lh,ev,pr,ust,vst)
    implicit none
    real(idx),allocatable :: sw(:,:),lw(:,:),sh(:,:),lh(:,:)
    real(idx),allocatable :: ev(:,:),pr(:,:)
    real(idx),allocatable :: ust(:,:),vst(:,:)
    deallocate(sw); deallocate(lw)
    deallocate(sh); deallocate(lh)
    deallocate(ev); deallocate(pr)
    deallocate(ust) ; deallocate(vst)
  end subroutine deallocate_atm_2d_arrays
  subroutine deallocate_mix_3d_arrays(bvf,shear,km,kt)
    implicit none
    real(idx),allocatable :: bvf(:,:,:),shear(:,:,:),km(:,:,:),kt(:,:,:)
    deallocate(bvf) ; deallocate(shear)
    deallocate(km) ;  deallocate(kt)
  end subroutine deallocate_mix_3d_arrays
  subroutine create_hist_rho(fname,nlon,nlat,nlev,lon,lat,lev, &
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
    write(*,*) "Create hist_rho=",trim(fname)    
  end subroutine create_hist_rho
  ! Output q-lev file
  subroutine create_hist_q(fname,nlon,nlat,nlev,ntime,lon,lat,z_q,time, &
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
    write(*,*) "Create hist_q=",trim(fname)    
  end subroutine create_hist_q
  subroutine output_hist_rho(fname,nlon,nlat,nlev,itime,&
       & out_sw,out_lw,out_sh,out_lh,out_ev,out_pr,&
       & out_ust,out_vst,out_temp,out_salt,out_u,out_v,&
       & sw,lw,sh,lh,ev,pr,ust,vst,&
       & temp,salt,u,v)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_sw,out_lw,out_sh,out_lh
    logical,intent(in) :: out_ev,out_pr,out_ust,out_vst
    logical,intent(in) :: out_temp,out_salt,out_u,out_v
    real(idx),intent(in) :: sw(nlon,nlat),lw(nlon,nlat),sh(nlon,nlat),lh(nlon,nlat)
    real(idx),intent(in) :: ev(nlon,nlat),pr(nlon,nlat),ust(nlon,nlat),vst(nlon,nlat)
    real(idx),intent(in) :: temp(nlon,nlat,nlev),salt(nlon,nlat,nlev),u(nlon,nlat,nlev),v(nlon,nlat,nlev)
    real(idx) :: tmp_3D(1:nlon,1:nlat,1)
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    tmp_3D = 0.0_idx
    tmp_4D = 0.0_idx
    if (out_sw .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=sw(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"sw",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_lw .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=lw(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"lw",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_sh .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=sh(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"sh",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_lh .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=lh(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"lh",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_ev .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=ev(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"ev",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_pr .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=pr(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"pr",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_ust .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=ust(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"ust",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_vst .eqv. .TRUE.) then
       tmp_3D(1:nlon,1:nlat,1)=vst(1:nlon,1:nlat)
       call writenet_wv(trim(fname),"vst",(/1,1,itime/),(/nlon,nlat,itime/),&
            & tmp_3D)
    end if
    if (out_temp .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=temp(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"temp",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_salt .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=salt(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"salt",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_u .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=u(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"u",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_v .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=v(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"v",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
  end subroutine output_hist_rho
  ! Output
  subroutine output_hist_q(fname,nlon,nlat,nlev,itime,&
       & out_bvf,out_shear,out_qq,out_l,out_km,out_kt,bvf,shear,&
#if defined(NNF) | defined(KC) | defined(bd)       
       & qq,l,&
#endif       
       & km,kt)
    implicit none
    character(len=*),intent(in) :: fname
    integer,intent(in) :: nlon,nlat,nlev,itime
    logical,intent(in) :: out_bvf,out_shear,out_qq,out_l
    logical,intent(in) :: out_km,out_kt
    real(idx),intent(in) :: bvf(nlon,nlat,nlev),shear(nlon,nlat,nlev)
#if defined(NNF) | defined(KC) | defined(bd)
    real(idx),intent(in) :: qq(nlon,nlat,nlev),l(nlon,nlat,nlev)
#endif
    real(idx),intent(in) :: km(nlon,nlat,nlev),kt(nlon,nlat,nlev)
    real(idx) :: tmp_4D(1:nlon,1:nlat,1:nlev,1)
    tmp_4D=0.0_idx
    if (out_bvf .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=bvf(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"bvf",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_shear .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=shear(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"shear",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
#if defined(NNF) | defined(KC) | defined(bd)
    if (out_qq .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=qq(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"qq",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_l .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=l(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"l",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
#endif
    if (out_km .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=km(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"km",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
    if (out_kt .eqv. .TRUE.) then
       tmp_4D(1:nlon,1:nlat,1:nlev,1)=kt(1:nlon,1:nlat,1:nlev)
       call writenet_wv(trim(fname),"kt",(/1,1,1,itime/),(/nlon,nlat,nlev,itime/),&
            & tmp_4D)
    end if
  end subroutine output_hist_q
end module mod_hist
