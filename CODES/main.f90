program main
  ! Main source program of mixed layer model
  use param
  use arrays_sub
  use input_files
  use mod_adv
  use mod_atm
  use mod_clm
  use mod_hist
  use mod_avg
  use mod_diags
  use mod_rst
  use mod_init
  use bulk_sub
  use eos_sub
  use solve_diag
  !$ use omp_lib
#include 'CPPLISTS.h'
#if defined(NNF)
  use nnf_sub
#endif
#if defined(KC)
  use kc_sub
#endif
#if defined(KPP)
  use kpp_sub
#endif
  implicit none
  character(maxlen) :: fname_grid
  integer :: nlon,nlat,nz,ntime
  ! Time step parameters
  real(idx) :: total_time,elapsed_time,tmp
  !===========================================
  ! Output setting
  integer :: iavg_count,itdiag_count,isdiag_count,iuvdiag_count
#if defined(NNF) | defined(KC) | defined(bd)
  integer :: iqqdiag_count
#endif
  integer :: ix,iy,iz,it
  !$ double precision st, en
  !===========================================
  namelist/input_grid/fname_grid
  !$ st = omp_get_wtime()
  !===========================================!
  ! Read namelists                            !
  !===========================================!
!  open(unit=nmlf_master,file=trim(namelist_master))
 ! read(unit=nmlf_master,nml=master)
  read(*,nml=master)
  !===========================================!
  ! Load parameters                           !
  !===========================================!
  open(unit=nmlf_set,file=trim(namelist_set))  
  read(unit=nmlf_set,nml=parameters)
  close(nmlf_set)
  tau_temp=tau_temp_day*day_to_sec; 
  tau_salt=tau_salt_day*day_to_sec; 
  tau_u=tau_u_day*day_to_sec; 
  tau_v=tau_v_day*day_to_sec; 
  !===========================================!
  ! Get information on time                   !
  !===========================================!
  open(unit=nmlf_io,file=trim(namelist_io))
  read(unit=nmlf_io,nml=time)
  call calendar_cal_length_ymdhms(start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int,-10000,tmp)
  total_time=int(tmp)
  ntime=int(total_time/dt)
  !===========================================!
  ! Read oceanic grid                         !
  !===========================================!
  read(unit=nmlf_io,nml=input_grid)
  write(*,*) "Grid=",trim(fname_grid)
  call read_grid_file(fname_grid,nlon,nlat,nz,lon_grd,lat_grd,cor,z_rho,z_q)
  call allocate_atm_2d_arrays(nlon,nlat,sw,lw,sh,lh,ev,pr,ust,vst)
  !===========================================!
  ! Read initial contidion                    !
  !===========================================!
  read(unit=nmlf_io,nml=input_init)
  call read_init_file(fname_init_temp,varname_init_temp,istart_temp,nlon,nlat,nz,temp)
  write(*,*) "Init temp =",trim(fname_init_temp)
  call read_init_file(fname_init_salt,varname_init_salt,istart_salt,nlon,nlat,nz,salt)
  write(*,*) "Init salt = ",trim(fname_init_salt)
  call read_init_file(fname_init_u,varname_init_u,istart_u,nlon,nlat,nz,u)
  write(*,*) "Init u = ",trim(fname_init_u)
  call read_init_file(fname_init_v,varname_init_v,istart_v,nlon,nlat,nz,v)
  write(*,*) "Init v = ",trim(fname_init_v)
  ! Initialize turbulent field-----------------------------------------------
#if defined(NNF) | defined(KC) | defined(bd)
    call read_init_file(fname_init_qq,varname_init_qq,istart_qq,nlon,nlat,nz,qq)
    write(*,*) "Init qq = ",trim(fname_init_qq)
    call read_init_file(fname_init_l,varname_init_l,istart_l,nlon,nlat,nz,l)
    write(*,*) "Init l =",trim(fname_init_l)
#endif
  !===========================================!
  ! Read atmospheric forcing file             !
  !===========================================!
    read(unit=nmlf_io,nml=input_atm_param)
    allocate(fnames_pr(nfile_pr));
    allocate(fnames_sw(nfile_sw)); allocate(fnames_lw(nfile_lw))
    allocate(fnames_ta(nfile_ta)); allocate(fnames_qa(nfile_qa))
    allocate(fnames_uw(nfile_uw)); allocate(fnames_vw(nfile_vw))
#if defined(WSPEED_IN)
    allocate(fnames_ws(nfile_ws))
#endif
    read(unit=nmlf_io,nml=input_atm_io)
    ! Precipitation
  call read_3D_files(fnames_pr,varname_pr,ntime_pr,time_pr,pr_in,start_yymmdd_int,start_hhmmss_int)
  ! Radiation
  write(*,*) trim(fnames_sw(1))
  call read_3D_files(fnames_sw,varname_sw,ntime_sw,time_sw,sw_in,start_yymmdd_int,start_hhmmss_int)
  call read_3D_files(fnames_lw,varname_lw,ntime_lw,time_lw,lw_in,start_yymmdd_int,start_hhmmss_int)
  call read_3D_files(fnames_ta,varname_ta,ntime_ta,time_ta,ta_in,start_yymmdd_int,start_hhmmss_int)
  call read_3D_files(fnames_qa,varname_qa,ntime_qa,time_qa,qa_in,start_yymmdd_int,start_hhmmss_int)
  call read_3D_files(fnames_uw,varname_uw,ntime_uw,time_uw,uw_in,start_yymmdd_int,start_hhmmss_int)
  call read_3D_files(fnames_vw,varname_vw,ntime_vw,time_vw,vw_in,start_yymmdd_int,start_hhmmss_int)
#if defined(WSPEED_IN)
  call read_3D_files(fnames_ws,varname_ws,ntime_ws,time_ws,ws_in,start_yymmdd_int,start_hhmmss_int)  
#endif
  write(*,*) "*******************************************************"
  write(*,*) " Finish reading atmospheric file"
  write(*,*) "*******************************************************"
  !===========================================!
  ! Read climatology file                     !
  !===========================================!
  read(unit=nmlf_io,nml=input_clm_param)
  allocate(fnames_clm_temp(nfile_clm_temp));  allocate(fnames_clm_salt(nfile_clm_salt))
  allocate(fnames_clm_u(nfile_clm_u));  allocate(fnames_clm_v(nfile_clm_v))
  read(unit=nmlf_io,nml=input_clm_io)
  if (use_clm_temp .eqv. .TRUE.) then
     call read_4D_files(fnames_clm_temp,varname_clm_temp,ntime_clm_temp,time_clm_temp,&
          & temp_clm,start_yymmdd_int,start_hhmmss_int)
     allocate(temp_clm_1d(nz))
  end if
  if (use_clm_salt .eqv. .TRUE.) then
     call read_4D_files(fnames_clm_salt,varname_clm_salt,ntime_clm_salt,time_clm_salt,&
          & salt_clm,start_yymmdd_int,start_hhmmss_int)
     allocate(salt_clm_1d(nz))
  end if
  if (use_clm_u .eqv. .TRUE.) then
     call read_4D_files(fnames_clm_u,varname_clm_u,ntime_clm_u,time_clm_u,&
          & u_clm,start_yymmdd_int,start_hhmmss_int)
     allocate(u_clm_1d(nz))
  end if
  if (use_clm_v .eqv. .TRUE.) then  
     call read_4D_files(fnames_clm_v,varname_clm_v,ntime_clm_v,time_clm_v,&
          & v_clm,start_yymmdd_int,start_hhmmss_int)
     allocate(v_clm_1d(nz))
  end if
  !===========================================!
  ! Read advection file                       !
  !===========================================!
  read(unit=nmlf_io,nml=input_adv_param)
  allocate(fnames_adv_temp(nfile_adv_temp));  allocate(fnames_adv_salt(nfile_adv_salt))
  allocate(fnames_adv_u(nfile_adv_u));  allocate(fnames_adv_v(nfile_adv_v))
  read(unit=nmlf_io,nml=input_adv_io)
  if (use_adv_temp .eqv. .TRUE.) then
     call read_4D_files(fnames_adv_temp,varname_adv_temp,ntime_adv_temp,time_adv_temp,&
          & temp_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  if (use_adv_salt .eqv. .TRUE.) then
     call read_4D_files(fnames_adv_salt,varname_adv_salt,ntime_adv_salt,time_adv_salt,&
          & salt_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  if (use_adv_u .eqv. .TRUE.) then
     call read_4D_files(fnames_adv_u,varname_adv_u,ntime_adv_u,time_adv_u,&
          & u_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  if (use_adv_v .eqv. .TRUE.) then  
     call read_4D_files(fnames_adv_v,varname_adv_v,ntime_adv_v,time_adv_v,&
          & v_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  !=================================================
  ! Allocate arrays necessary for calculation
  !==================================================
  allocate(bvf(nlon,nlat,nz)) ; allocate(shear(nlon,nlat,nz))
  allocate(km(nlon,nlat,nz)) ;  allocate(kt(nlon,nlat,nz))
  allocate(temp_next(1:nz)) ; allocate(salt_next(1:nz))
  allocate(u_next(1:nz)) ; allocate(v_next(1:nz))
  temp_next=0.0_idx; salt_next=0.0_idx; u_next=0.0_idx; v_next=0.0_idx
#if defined(NNF) | defined(KC) | defined(bd)
  allocate(qq_next(1:nz)); allocate(l_next(1:nz))
  qq_next=0.0_idx ; l_next=0.0_idx
#endif
  allocate(bvf_1d(1:(nz-1))); allocate(shear_1d(1:(nz-1)))
  allocate(km_1d(1:(nz-1))); allocate(kt_1d(1:(nz-1))) ; allocate(ks_1d(1:(nz-1)))
  bvf_1d=0.0_idx; shear_1d=0.0_idx ; km_1d=0.0_idx; kt_1d=0.0_idx; ks_1d=0.0_idx
#if defined(NNF) | defined(KC) | defined(bd)
  allocate(kq_1d(1:(nz-1))); kq_1d=0.0_idx
# endif
  allocate(penet(nz))
  allocate(temp_adv_1d(nz)) ; allocate(salt_adv_1d(nz))
  allocate(u_adv_1d(nz)) ;   allocate(v_adv_1d(nz))
  temp_next=0.0_idx ; salt_next =0.0_idx;  u_next = 0.0_idx  ; v_next = 0.0_idx
#if defined(NNF) | defined(KC) | defined(bd)
  qq_next = 0.0_idx ; l_next=0.0_idx
#endif
  penet=0.0_idx
  !======================!
  ! Prepare output files !
  !======================!
  !*********************************
  ! History (snapshot) file
  !*********************************
  read(unit=nmlf_io,nml=output_hist_rho_flag)
  call create_hist_rho(fname_out_hist_rho,nlon,nlat,nz,lon_grd,lat_grd,z_rho,&
       & start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int, &
       & out_hist_flag,out_hist_int,&
       & hout_temp,hout_salt,hout_u,hout_v,&
       & hout_sw,hout_lw,hout_sh,hout_lh,&
       & hout_ev,hout_pr,hout_ust,hout_vst,&
       & missing_value,istep_hist,ntime_hist,time_hist)
  ! Q-file)
  read(unit=nmlf_io,nml=output_hist_q_flag)
  call create_hist_q(fname_out_hist_q,nlon,nlat,nz,ntime_hist,lon_grd,lat_grd,z_q,&
       & time_hist, start_yymmdd_int,start_hhmmss_int,out_hist_flag, &
       & hout_bvf,hout_shear,hout_qq,hout_l,hout_km,hout_kt,missing_value)
  !*********************************!
  ! Average file                    !
  !*********************************!
  allocate(temp_avg(nlon,nlat,nz)) ; temp_avg=0.0_idx
  allocate(salt_avg(nlon,nlat,nz)) ; salt_avg=0.0_idx
  allocate(u_avg(nlon,nlat,nz)) ; u_avg=0.0_idx
  allocate(v_avg(nlon,nlat,nz)) ; v_avg=0.0_idx
#if defined(NNF) | defined(KC) | defined(bd)
  allocate(qq_avg(nlon,nlat,nz)) ; qq_avg=0.0_idx
  allocate(l_avg(nlon,nlat,nz)) ; l_avg=0.0_idx
#endif  
  allocate(bvf_avg(nlon,nlat,nz)) ; bvf_avg=0.0_idx
  allocate(shear_avg(nlon,nlat,nz)); shear_avg=0.0_idx
  allocate(km_avg(nlon,nlat,nz)) ; km_avg=0.0_idx
  allocate(kt_avg(nlon,nlat,nz)); kt_avg=0.0_idx
  call allocate_atm_2d_arrays(nlon,nlat,sw_avg,lw_avg,sh_avg,lh_avg,ev_avg,pr_avg,ust_avg,vst_avg)
  read(unit=nmlf_io,nml=output_avg_rho_flag)
  call create_avg_rho(fname_out_avg_rho,nlon,nlat,nz,lon_grd,lat_grd,z_rho,&
       & start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int, &
       & out_avg_flag,out_avg_int,&
       & hout_temp,hout_salt,hout_u,hout_v,&
       & hout_sw,hout_lw,hout_sh,hout_lh,&
       & hout_ev,hout_pr,hout_ust,hout_vst,&
       & missing_value,istep_avg,ntime_avg,time_avg)
  read(unit=nmlf_io,nml=output_avg_q_flag)
  ! File  (q-file)
  call create_avg_q(fname_out_avg_q,nlon,nlat,nz,ntime_avg,lon_grd,lat_grd,z_q,&
       & time_avg, start_yymmdd_int,start_hhmmss_int,out_avg_flag,& 
       & aout_bvf,aout_shear,aout_qq,aout_l,aout_km,aout_kt,missing_value)
  !*********************************
  ! Diagnostic (temperature) file
  !*********************************
  read(unit=nmlf_io,nml=output_temp_diag_flag)
  call allocate_tdiag(nlon,nlat,nz,temp_rate,temp_vdiff,temp_penet,temp_relax)
  call create_tdiag(fname_out_tdiag,nlon,nlat,nz,lon_grd,lat_grd,z_rho,&
       & start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int, &
       & out_tdia_flag,out_tdia_int,dout_temp_rate,dout_temp_vdiff,dout_temp_penet,dout_temp_relax,&
       & missing_value,istep_tdia,ntime_tdia,time_tdia)
  write(*,*) "Finish preparing "//trim(fname_out_tdiag)
  !*********************************
  ! Diagnostic (salinity) file
  !*********************************
  read(unit=nmlf_io,nml=output_salt_diag_flag)
  call allocate_sdiag(nlon,nlat,nz,salt_rate,salt_vdiff,salt_relax)
  call create_sdiag(fname_out_sdiag,nlon,nlat,nz,lon_grd,lat_grd,z_rho,&
       & start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int, &
       & out_sdia_flag,out_sdia_int,dout_salt_rate,dout_salt_vdiff,dout_salt_relax,&
       & missing_value,istep_sdia,ntime_sdia,time_sdia)
  write(*,*) "Finish preparing "//trim(fname_out_sdiag)
  !*********************************
  ! Diagnostic (uv) file
  !*********************************
  read(unit=nmlf_io,nml=output_uv_diag_flag)
  call allocate_uvdiag(nlon,nlat,nz,u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
  call create_uvdiag(fname_out_uvdiag,nlon,nlat,nz,lon_grd,lat_grd,z_rho,&
       & start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int, &
       & out_sdia_flag,out_sdia_int,dout_uv_rate,dout_uv_cor,dout_uv_vdiff,dout_uv_relax,&
       & missing_value,istep_uvdia,ntime_uvdia,time_uvdia)
  !*********************************
  ! Diagnostic (qq) file
  !*********************************
  read(unit=nmlf_io,nml=output_qq_diag_flag)
#if defined(NNF) | defined(KC) | defined(bd)
  call allocate_qqdiag(nlon,nlat,nz,qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp)
  call create_qqdiag(fname_out_qqdiag,nlon,nlat,nz,lon_grd,lat_grd,z_q,&
       & start_yymmdd_int,start_hhmmss_int,end_yymmdd_int,end_hhmmss_int, &
       & out_qqdia_flag,out_qqdia_int,dout_qq_rate,dout_qq_vdiff,dout_qq_sp,&
       & dout_qq_bp,dout_qq_disp,&
       & missing_value,istep_qqdia,ntime_qqdia,time_qqdia)
#endif  
  close(nmlf_io)

  !==========================!
  !  Initialize output flag  !
  !==========================!
  ihist = 1 ; iavg = 1 ; itdia = 1 ; isdia = 1 ; iuvdia = 1;
  iavg_count=0 ; itdiag_count=0 ; isdiag_count=0; iuvdiag_count = 0
#if defined(NNF) | defined(KC) | defined(bd)
  iqqdia=1
  iqqdiag_count=0
#endif
  !==========================!
  !  Start point of loop     !
  !==========================!
  write(*,*) "nlon,nlat,nz=",nlon,nlat,nz
  write(*,*) "Total timestep=",ntime
!  ntime=1
  do it = 1,ntime     
     iavg_count = iavg_count + 1
     itdiag_count = itdiag_count + 1
     isdiag_count = isdiag_count + 1
     iuvdiag_count = iuvdiag_count + 1
#if defined(NNF) | defined(KC) | defined(bd)
     iqqdiag_count = iqqdiag_count + 1
#endif
     elapsed_time=dt*it
     ! Time interpolation
     call time_wgt2(time_pr,elapsed_time,ind1_pr,ind2_pr,wgt1_pr,wgt2_pr)
     call time_wgt2(time_sw,elapsed_time,ind1_sw,ind2_sw,wgt1_sw,wgt2_sw)
     call time_wgt2(time_lw,elapsed_time,ind1_lw,ind2_lw,wgt1_lw,wgt2_lw)
     call time_wgt2(time_ta,elapsed_time,ind1_ta,ind2_ta,wgt1_ta,wgt2_ta)
     call time_wgt2(time_qa,elapsed_time,ind1_qa,ind2_qa,wgt1_qa,wgt2_qa)
     call time_wgt2(time_uw,elapsed_time,ind1_uw,ind2_uw,wgt1_uw,wgt2_uw)
     call time_wgt2(time_vw,elapsed_time,ind1_vw,ind2_vw,wgt1_vw,wgt2_vw)
#if defined(WSPEED_IN)
     call time_wgt2(time_ws,elapsed_time,ind1_ws,ind2_ws,wgt1_ws,wgt2_ws)
#endif
     ! temp
     if (use_clm_temp .eqv. .TRUE.) then
        call time_wgt(time_clm_temp,elapsed_time,ind1_clm_temp,ind2_clm_temp,wgt1_clm_temp,wgt2_clm_temp)
     end if
     ! salt
     if (use_clm_salt .eqv. .TRUE.) then
        call time_wgt(time_clm_salt,elapsed_time,ind1_clm_salt,ind2_clm_salt,wgt1_clm_salt,wgt2_clm_salt)
     end if
     ! u
     if (use_clm_u .eqv. .TRUE.) then
        call time_wgt(time_clm_u,elapsed_time,ind1_clm_u,ind2_clm_u,wgt1_clm_u,wgt2_clm_u)
     end if
     if (use_clm_v .eqv. .TRUE.) then
        call time_wgt(time_clm_v,elapsed_time,ind1_clm_v,ind2_clm_v,wgt1_clm_v,wgt2_clm_v)
     end if
     ! temp
     if (use_adv_temp .eqv. .TRUE.) then        
        call time_wgt2(time_adv_temp,elapsed_time,ind1_adv_temp,ind2_adv_temp,wgt1_adv_temp,wgt2_adv_temp)
     end if
     ! salt
     if (use_adv_salt .eqv. .TRUE.) then
        call time_wgt2(time_adv_salt,elapsed_time,ind1_adv_salt,ind2_adv_salt,wgt1_adv_salt,wgt2_adv_salt)
     end if
     ! u
     if (use_adv_u .eqv. .TRUE.) then
        call time_wgt2(time_adv_u,elapsed_time,ind1_adv_u,ind2_adv_u,wgt1_adv_u,wgt2_adv_u)
     end if
     if (use_adv_v .eqv. .TRUE.) then
!        call time_wgt(time_adv_v,elapsed_time,ind1_adv_v,ind2_adv_v,wgt1_adv_v,wgt2_adv_v)
        call time_wgt2(time_adv_v,elapsed_time,ind1_adv_v,ind2_adv_v,wgt1_adv_v,wgt2_adv_v)

     end if

     ! Lon-lat loop
     !$omp parallel do private(ld,ta,qa,uwind,vwind,ws,heat_solar,heat_no_solar,e_p,ssflux,&
     !$omp & penet,bvf_1d,shear_1d,km_1d,kt_1d,ks_1d,kq_1d,&
     !$omp & temp_adv_1d,salt_adv_1d,u_adv_1d,v_adv_1d,&
     !$omp & temp_clm_1d,salt_clm_1d,u_clm_1d,v_clm_1d,&
     !$omp & temp_next,salt_next,u_next,v_next,qq_next,l_next)
     do iy = 1,nlat
        do ix = 1,nlon
           !=============================!
           ! Specify boundary conditions !
           !=============================!
           if (abs(temp(ix,iy,nz)) .lt. missing_value) then              
              ! Precipitation
              pr(ix,iy)=set_data(ind1_pr,ind2_pr,wgt1_pr,wgt2_pr,pr_in(ix,iy,:))
              ! Radiation
              sw(ix,iy)=set_data(ind1_sw,ind2_sw,wgt1_sw,wgt2_sw,sw_in(ix,iy,:))
# if defined(LONGWAVE_DOWN)
              ld=set_data(ind1_lw,ind2_lw,wgt1_lw,wgt2_lw,lw_in(ix,iy,:))
              lw(ix,iy)=ld_to_lw(ld,temp(ix,iy,nz))
# else
              lw(ix,iy)=set_data(ind1_lw,ind2_lw,wgt1_lw,wgt2_lw,lw_in(ix,iy,:))
# endif
              ! Wind stress
# if defined(WSTRESS_IN)
              ust(ix,iy)=set_data(ind1_uw,ind2_uw,wgt1_uw,wgt2_uw,uw_in(ix,iy,:))
              vst(ix,iy)=set_data(ind1_vw,ind2_vw,wgt1_vw,wgt2_vw,vw_in(ix,iy,:))
# else
              uwind=set_data(ind1_uw,ind2_uw,wgt1_uw,wgt2_uw,uw_in(ix,iy,:))
              vwind=set_data(ind1_vw,ind2_vw,wgt1_vw,wgt2_vw,vw_in(ix,iy,:))
#  if defined(WSPEED_IN)
              ws=set_data(ind1_ws,ind2_ws,wgt1_ws,wgt2_ws,ws_in(ix,iy,:))
#  else
              ws=sqrt(uwind*uwind+vwind*vwind)
#  endif
              ! SPEED
              call cal_ws(uwind,vwind,ws,ust(ix,iy),vst(ix,iy))
# endif
              ! Turbulent fluxes
#if defined(BULK_FLUXES)
              ta=set_data(ind1_ta,ind2_ta,wgt1_ta,wgt2_ta,ta_in(ix,iy,:))
              qa=set_data(ind1_qa,ind2_qa,wgt1_qa,wgt2_qa,qa_in(ix,iy,:))
#  if defined(WSPEED_IN)
              ws=set_data(ind1_ws,ind2_ws,wgt1_ws,wgt2_ws,ws_in(ix,iy,:))
#  else
              ws=sqrt(uwind*uwind+vwind*vwind)
#  endif
# if defined (BULK_KARA05)
              call bulk_kara05(ta,qa,temp(ix,iy,nz),ws,sh(ix,iy),lh(ix,iy),ev(ix,iy))
# else 
              call bulk_ncep(ta,qa,temp(ix,iy,nz),ws,sh(ix,iy),lh(ix,iy),ev(ix,iy))
# endif
#else
              sh(ix,iy)=set_data(ind1_ta,ind2_ta,wgt1_ta,wgt2_ta,ta_in(ix,iy,:))
              lh(ix,iy)=set_data(ind1_qa,ind2_qa,wgt1_qa,wgt2_qa,qa_in(ix,iy,:))
              ev(ix,iy)= lh(ix,iy)* lh_to_ev
#endif
              !====================!
              ! Set shortwave flux !                     
              !====================!
#if defined(DIURNAL)
              heat_solar = diurnal_sw(sw(ix,iy),elapsed_time)
#else
              heat_solar = sw(ix,iy)
#endif
              call heat_absorb(nz,z_q,heat_solar,beta1,beta2,r_long,penet)
              !=========================!
              ! Set non-solar heat flux !
              !=========================!
              heat_no_solar = lw(ix,iy) + sh(ix,iy) + lh(ix,iy)
              !=====================!
              ! Set freshwater flux !
              !=====================!
              e_p = abs(ev(ix,iy)) - abs(pr(ix,iy))  ! evaporation minus precipitation
              ssflux=e_p * 35.0_idx

              !ssflux=e_p * salt(ix,iy,nz)
              !========================================!
              ! Calculate buoyancy frequency and shear !
              !========================================!
              bvf_1d(1:nz-1)=cal_bv(nz,z_rho,temp(ix,iy,1:nz),salt(ix,iy,1:nz))
              shear_1d(1:nz-1)=cal_shear(nz,z_rho(1:nz),u(ix,iy,1:nz),v(ix,iy,1:nz))
              !=========================================!
              ! Compute vertical viscocity coefficients !
              !=========================================!
#if defined(KC)
              call my25_vdiff(nz,bvf_1d,shear_1d,qq(ix,iy,:),l(ix,iy,:),&
                   & km_1d,kt_1d,ks_1d)
#endif
#if defined(NNF)
              call mynnf25_vdiff(nz,bvf_1d,shear_1d,qq(ix,iy,:),l(ix,iy,:),&
                   & km_1d,kt_1d,ks_1d)
#endif
#if defined(bd)
              call bd_vdiff(nz,bvf_1d,shear_1,qq(ix,iy,:),l(ix,iy,:),&
                   & km_1,kt_1d,ks_1d)
#endif
#if defined(KPP)
              call kpp_vdiff(nz,z_rho,z_q,temp(ix,iy,:),salt(ix,iy,:),u(ix,iy,:),v(ix,iy,:),&
                   & bvf_1d,shear_1d,km_1d,kt_1d,ks_1d,&
                   & ust(ix,iy),vst(ix,iy),heat_solar,heat_no_solar,e_p,beta1,beta2,r_long,cor(iy))
#endif
              !===========================================!
              ! Calculate turbulent kinetic energy(QQ)    !
              !===========================================!
#if defined(NNF) | defined(KC) | defined(bd)
              call cal_kq(nz,km_1d,kq_1d)
              call cal_next_qq_bud(dt,nz,z_rho,z_q,qq(ix,iy,:),bvf_1d(1:nz-1),shear_1d(1:nz-1) &
                   & ,l(ix,iy,:),kq_1d,km_1d,kt_1d,ust(ix,iy),vst(ix,iy),B_1,qq_next,&
                   & qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp)

#endif              
              !===========================================!
              ! Calculate turbulence length scale (L)     !
              !===========================================!
#if defined(NNF)
              call diag_l_mynnf(nz,z_q,bvf_1d,qq(ix,iy,:),l_next,&
                   & ust(ix,iy),vst(ix,iy),heat_solar,heat_no_solar,e_p,beta1,beta2,r_long)
#endif
#if defined(bd)
              call diag_l_bd(nz,z_q,bvf_1d,qq(ix,iy,:),l)
#endif              
#if defined(KC)
              call cal_next_l(dt,nz,z_rho,z_q,qq(ix,iy,:),bvf_1d(1:nz-1),shear_1d(1:nz-1) &
                   & ,l(ix,iy,:),kq_1d,km_1d,kt_1d,qq_next,B_1,E_1,E_2,l_next)
#endif
              !===========================================!
              ! Add background diffusion                  !
              !===========================================!
              do iz=1,nz-1
                 km_1d(iz) = km_1d(iz)  + nu ; kt_1d(iz) = kt_1d(iz)  + nu_t
                 ks_1d(iz) = ks_1d(iz)  + nu_s
              end do
              !=============================================!
              ! Calculate next values using implicit method !
              !=============================================!
              ! ** Temperature equation**
              if (use_adv_temp .eqv. .TRUE.) then
                 do iz = 1,nz
                    temp_adv_1d(iz)=set_data(ind1_adv_temp,ind2_adv_temp,wgt1_adv_temp,wgt2_adv_temp,&
                         & temp_adv(ix,iy,iz,:))
                 end do
              else                
                 temp_adv_1d = 0.0_idx
              end if
              if (use_clm_temp .eqv. .TRUE.) then
                 do iz = 1,nz
                    temp_clm_1d(iz)=set_data(ind1_clm_temp,ind2_clm_temp,wgt1_clm_temp,wgt2_clm_temp,&
                         & temp_clm(ix,iy,iz,:))
                    temp_adv_1d(iz)=temp_adv_1d(iz)+(temp_clm_1d(iz)-temp(ix,iy,iz))/tau_temp                                        
                 end do
              end if
              call cal_next_theta_bud(dt,nz,z_rho,z_q,temp(ix,iy,1:nz),kt_1d,penet,&
                   & heat_no_solar,temp_adv_1d,temp_next,&
                   & temp_rate(ix,iy,1:nz),temp_vdiff(ix,iy,1:nz),&
                   & temp_penet(ix,iy,1:nz),temp_relax(ix,iy,1:nz))
                 
              ! Salinity equation
              if (use_adv_salt .eqv. .TRUE.) then
                 do iz = 1,nz
                    salt_adv_1d(iz)=set_data(ind1_adv_salt,ind2_adv_salt,wgt1_adv_salt,wgt2_adv_salt,&
                         & salt_adv(ix,iy,iz,:))
                 end do
              else
                 salt_adv_1d = 0.0_idx                
              end if
              if (use_clm_salt .eqv. .TRUE.) then
                 do iz = 1,nz
                    salt_clm_1d(iz)=set_data(ind1_clm_salt,ind2_clm_salt,wgt1_clm_salt,wgt2_clm_salt,&
                         & salt_clm(ix,iy,iz,:))
                    salt_adv_1d(iz)=salt_adv_1d(iz)+(salt_clm_1d(iz)-salt(ix,iy,iz))/tau_salt
                 end do
              end if
              call cal_next_salt_bud(dt,nz,z_rho,z_q,salt(ix,iy,1:nz),&
                   & ks_1d,ssflux,salt_adv_1d,salt_next,&
                   & salt_rate(ix,iy,1:nz),salt_vdiff(ix,iy,1:nz),salt_relax(ix,iy,1:nz))

              ! Velocity equation
              if (use_adv_u .eqv. .TRUE.) then
                 do iz = 1,nz
                    u_adv_1d(iz)=set_data(ind1_adv_u,ind2_adv_u,wgt1_adv_u,wgt2_adv_u,&
                         & u_adv(ix,iy,iz,:))
                 end do
              else
                 u_adv_1d = 0.0_idx                
              end if
              if (use_adv_v .eqv. .TRUE.) then
                 do iz = 1,nz
                    v_adv_1d(iz)=set_data(ind1_adv_v,ind2_adv_v,wgt1_adv_v,wgt2_adv_v,&
                         & v_adv(ix,iy,iz,:))
                 end do
              else
                 v_adv_1d = 0.0_idx                
              end if
              if (use_clm_u .eqv. .TRUE.) then
                 do iz = 1,nz
                    u_clm_1d(iz)=set_data(ind1_clm_u,ind2_clm_u,wgt1_clm_u,wgt2_clm_u,&
                         & u_clm(ix,iy,iz,:))
                    u_adv_1d(iz)=u_adv_1d(iz)+(u_clm_1d(iz)-u(ix,iy,iz))/tau_u
                 end do
              end if
              if (use_clm_v .eqv. .TRUE.) then
                 do iz = 1,nz
                    v_clm_1d(iz)=set_data(ind1_clm_v,ind2_clm_v,wgt1_clm_v,wgt2_clm_v,&
                         & v_clm(ix,iy,iz,:))
                    v_adv_1d(iz)=v_adv_1d(iz)+(v_clm_1d(iz)-v(ix,iy,iz))/tau_v
                 end do
              end if
              call cal_next_uv_bud(dt,nz,z_rho,z_q,u(ix,iy,1:nz),v(ix,iy,1:nz),&
                   & km_1d,cor(iy),ust(ix,iy),vst(ix,iy),u_adv_1d,v_adv_1d,u_next,v_next,&
                   & u_rate(ix,iy,1:nz),u_cor(ix,iy,1:nz),u_vdiff(ix,iy,1:nz),u_relax(ix,iy,1:nz),&
                   & v_rate(ix,iy,1:nz),v_cor(ix,iy,1:nz),v_vdiff(ix,iy,1:nz),v_relax(ix,iy,1:nz))
              !================!
              ! Update values  !
              !================!
              do iz =1,nz
                 temp(ix,iy,iz)=temp_next(iz) ;  salt(ix,iy,iz)=salt_next(iz)
                 u(ix,iy,iz)=u_next(iz) ;  v(ix,iy,iz)=v_next(iz)
#if defined(NNF) | defined(KC) | defined(bd)
                 qq(ix,iy,iz) = qq_next(iz) ; l(ix,iy,iz) = l_next(iz)
#endif
              end do
              !===========!
              ! Averaging !
              !===========!
              sw_avg(ix,iy)=sw_avg(ix,iy)+sw(ix,iy)
              lw_avg(ix,iy)=lw_avg(ix,iy)+lw(ix,iy)
              sh_avg(ix,iy)=sh_avg(ix,iy)+sh(ix,iy)
              lh_avg(ix,iy)=lh_avg(ix,iy)+lh(ix,iy)
              ev_avg(ix,iy)=ev_avg(ix,iy)+ev(ix,iy)
              pr_avg(ix,iy)=pr_avg(ix,iy)+pr(ix,iy)
              ust_avg(ix,iy)=ust_avg(ix,iy)+ust(ix,iy)
              vst_avg(ix,iy)=vst_avg(ix,iy)+vst(ix,iy)
              temp_avg(ix,iy,1:nz)=temp_avg(ix,iy,1:nz)+temp(ix,iy,1:nz)
              salt_avg(ix,iy,1:nz)=salt_avg(ix,iy,1:nz)+salt(ix,iy,1:nz)
              u_avg(ix,iy,1:nz)=u_avg(ix,iy,1:nz)+u(ix,iy,1:nz)
              v_avg(ix,iy,1:nz)=v_avg(ix,iy,1:nz)+v(ix,iy,1:nz)
              bvf_avg(ix,iy,1:nz-1)=bvf_avg(ix,iy,1:nz-1)+bvf_1d(1:nz-1)
              shear_avg(ix,iy,1:nz-1)=shear_avg(ix,iy,1:nz-1)+shear_1d(1:nz-1)
#if defined(NNF) | defined(KC) | defined(bd)
              qq_avg(ix,iy,1:nz)=qq_avg(ix,iy,1:nz)+qq(ix,iy,1:nz)
              l_avg(ix,iy,1:nz)=l_avg(ix,iy,1:nz)+l(ix,iy,1:nz)
#endif
              km_avg(ix,iy,1:nz-1)=km_avg(ix,iy,1:nz-1)+km_1d(1:nz-1)
              kt_avg(ix,iy,1:nz-1)=kt_avg(ix,iy,1:nz-1)+kt_1d(1:nz-1)
           else
              temp(ix,iy,:) = missing_value ; salt(ix,iy,:) = missing_value
              u(ix,iy,:) = missing_value ; v(ix,iy,:) = missing_value    
              sw_avg(ix,iy)=sw_avg(ix,iy)+missing_value
              lw_avg(ix,iy)=lw_avg(ix,iy)+missing_value
              sh_avg(ix,iy)=sh_avg(ix,iy)+missing_value
              lh_avg(ix,iy)=lh_avg(ix,iy)+missing_value
              ev_avg(ix,iy)=ev_avg(ix,iy)+missing_value
              pr_avg(ix,iy)=pr_avg(ix,iy)+missing_value
              ust_avg(ix,iy)=ust_avg(ix,iy)+missing_value
              vst_avg(ix,iy)=vst_avg(ix,iy)+missing_value
              temp_avg(ix,iy,1:nz)=temp_avg(ix,iy,1:nz)+missing_value
              salt_avg(ix,iy,1:nz)=salt_avg(ix,iy,1:nz)+missing_value
              u_avg(ix,iy,1:nz)=u_avg(ix,iy,1:nz)+missing_value
              v_avg(ix,iy,1:nz)=v_avg(ix,iy,1:nz)+missing_value
#if defined(NNF) | defined(KC) | defined(bd)
              qq_avg(ix,iy,1:nz)=qq_avg(ix,iy,1:nz)+missing_value
              l_avg(ix,iy,1:nz)=l_avg(ix,iy,1:nz)+missing_value
#endif
              bvf_avg(ix,iy,1:nz-1)=bvf_avg(ix,iy,1:nz-1)+missing_value
              shear_avg(ix,iy,1:nz-1)=shear_avg(ix,iy,1:nz-1)+missing_value
              km_avg(ix,iy,1:nz-1)=km_avg(ix,iy,1:nz-1)+missing_value
              kt_avg(ix,iy,1:nz-1)=kt_avg(ix,iy,1:nz-1)+missing_value
              bvf(ix,iy,1:nz-1)=missing_value
              shear(ix,iy,1:nz-1)=missing_value
              km(ix,iy,1:nz-1)=missing_value
              kt(ix,iy,1:nz-1)=missing_value
           end if
           if (it .eq. istep_hist(ihist)) then
              bvf(ix,iy,1:nz-1) = bvf_1d(1:nz-1) ; shear(ix,iy,1:nz-1)=shear_1d(1:nz-1)
              km(ix,iy,1:nz-1) = km_1d(1:nz-1)   ; kt(ix,iy,1:nz-1) = kt_1d(1:nz-1)
           end if
        end do
     end do

     !==================!
     ! output procedure !
     !==================!
     if (it .eq. istep_hist(ihist)) then
        write(*,*) elapsed_time*sec_to_day, "Km maxlamda",0.5_idx*maxval(km)*dt/((z_q(nz)-z_q(nz-1))**2)
        call output_hist_rho(fname_out_hist_rho,nlon,nlat,nz,ihist, &
             & hout_sw,hout_lw,hout_sh,hout_lh,&
             & hout_ev,hout_pr,hout_ust,hout_vst,&
             & hout_temp,hout_salt,hout_u,hout_v,&
             & sw,lw,sh,lh,ev,pr,ust,vst,&
             & temp,salt,u,v)
        call output_hist_q(fname_out_hist_q,nlon,nlat,nz,ihist &
             & ,hout_bvf,hout_shear,hout_qq,hout_l,hout_km,hout_kt&
             & ,bvf,shear,&
#if defined(NNF) | defined(KC) | defined(bd)             
             & qq,l,&
#endif
             & km,kt)
        ihist=ihist+1
        ihist=min(ihist,ntime_hist)
     end if
     ! Average file
     if (it .eq. istep_avg(iavg)) then
        call output_avg_rho(fname_out_avg_rho,nlon,nlat,nz,iavg, &
             & aout_sw,aout_lw,aout_sh,aout_lh,&
             & aout_ev,aout_pr,aout_ust,aout_vst,&
             & aout_temp,aout_salt,aout_u,aout_v,&
             & sw_avg,lw_avg,sh_avg,lh_avg,ev_avg,pr_avg,ust_avg,vst_avg,&
             & temp_avg,salt_avg,u_avg,v_avg,iavg_count)
        call output_avg_q(fname_out_avg_q,nlon,nlat,nz,iavg &
             & ,aout_bvf,aout_shear,aout_qq,aout_l,aout_km,aout_kt&
             & ,bvf_avg,shear_avg &
#if defined(NNF) | defined(KC) | defined(bd)                          
             & ,qq_avg,l_avg &
#endif
             & ,km_avg,kt_avg,iavg_count)
        iavg=iavg+1
        iavg=min(iavg,ntime_avg)
        iavg_count=0
     end if
     ! Diagnostic file for temperature
     if (it .eq. istep_tdia(itdia)) then
        call output_diag_temp(fname_out_tdiag,nlon,nlat,nz,itdia &
             & ,dout_temp_rate,dout_temp_vdiff,dout_temp_penet,dout_temp_relax&
             & ,temp_rate,temp_vdiff,temp_penet,temp_relax,itdiag_count)
        itdia=itdia+1
        itdia=min(itdia,ntime_tdia)
        itdiag_count=0
     end if
     ! Diagnostic file for salinity
     if (it .eq. istep_sdia(isdia)) then
        call output_diag_salt(fname_out_sdiag,nlon,nlat,nz,isdia &
             & ,dout_salt_rate,dout_salt_vdiff,dout_salt_relax&
             & ,salt_rate,salt_vdiff,salt_relax,isdiag_count)
        isdia=isdia+1
        isdia=min(isdia,ntime_sdia)
        isdiag_count=0
     end if
     ! Diagnostic file for velocity
     if (it .eq. istep_uvdia(iuvdia)) then
        call output_diag_uv(fname_out_uvdiag,nlon,nlat,nz,iuvdia, &
             & dout_uv_rate,dout_uv_cor,dout_uv_vdiff,dout_uv_relax,&
             & u_rate,u_cor,u_vdiff,u_relax,&
             & v_rate,v_cor,v_vdiff,v_relax,iuvdiag_count)
        iuvdia=iuvdia+1
        iuvdia=min(iuvdia,ntime_uvdia)
        iuvdiag_count=0
     end if
#if defined(NNF) | defined(KC) | defined(bd)
     ! Diagnostic file for TKE
     if (it .eq. istep_qqdia(iqqdia)) then
        call output_diag_qq(fname_out_qqdiag,nlon,nlat,nz,iqqdia &
             & ,dout_qq_rate,dout_qq_vdiff,dout_qq_sp,dout_qq_bp,dout_qq_disp&
             & ,qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp,iqqdiag_count)
        iqqdia=iqqdia+1
        iqqdia=min(iqqdia,ntime_qqdia)
        iqqdiag_count=0
     end if
#endif
  end do
  it=ntime
  ! Output restart file
  open(unit=nmlf_io,file=trim(namelist_io))
  read(unit=nmlf_io,nml=output_rst_flag)
  close(nmlf_io)
  !$ en = omp_get_wtime()
  !$ write(*,*) "Elapsed time :", en-st
  call prepare_output_rst_rho(fname_out_rst_rho,nlon,nlat,nz,lon_grd,lat_grd,z_rho,&
       & it*dt/(60.0*60.0*24.0), start_yymmdd_int,start_hhmmss_int,1,missing_value)
  call output_rst_rho(fname_out_rst_rho,nlon,nlat,nz,temp,salt,u,v)
  call prepare_output_rst_q(fname_out_rst_q,nlon,nlat,nz,lon_grd,lat_grd,z_q,&
       & it*dt/(60.0*60.0*24.0), start_yymmdd_int,start_hhmmss_int,1,missing_value)
#if defined(NNF) | defined(KC) | defined(bd)
  call output_rst_q(fname_out_rst_q,nlon,nlat,nz,qq,l)
#endif
  ! Deallocation
  deallocate(lon_grd); deallocate(lat_grd); deallocate(z_rho)
  deallocate(z_q)
  deallocate(temp);deallocate(salt); deallocate(u); deallocate(v)
  deallocate(temp_next);deallocate(salt_next); deallocate(u_next); deallocate(v_next)
  deallocate(temp_avg);deallocate(salt_avg); deallocate(u_avg); deallocate(v_avg)
#if defined(NNF) | defined(KC) | defined(bd)
  deallocate(qq); deallocate(l)
  deallocate(qq_next); deallocate(l_next)
  deallocate(qq_avg); deallocate(l_avg)
  deallocate(kq_1d)
#endif
  deallocate(bvf);deallocate(shear); deallocate(km) ;deallocate(kt)
  deallocate(bvf_1d);deallocate(shear_1d); deallocate(km_1d)
  deallocate(kt_1d); deallocate(ks_1d)
  call deallocate_atm_2d_arrays(sw,lw,sh,lh,ev,pr,ust,vst)
  call deallocate_tdiag(temp_rate,temp_vdiff,temp_penet,temp_relax)
  deallocate(time_tdia); deallocate(istep_tdia)
  call deallocate_sdiag(salt_rate,salt_vdiff,salt_relax)
  deallocate(time_sdia); deallocate(istep_sdia)
  call deallocate_uvdiag(u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
  deallocate(time_uvdia); deallocate(istep_uvdia)
#if defined(NNF) | defined(KC) | defined(bd)
  call deallocate_qqdiag(qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp)
  deallocate(time_qqdia); deallocate(istep_qqdia)
#endif
  ! Atmospheric data
  deallocate(time_pr); deallocate(pr_in); deallocate(fnames_pr)
  deallocate(time_sw); deallocate(sw_in); deallocate(fnames_sw)
  deallocate(time_lw); deallocate(lw_in); deallocate(fnames_lw)
  deallocate(time_ta); deallocate(ta_in); deallocate(fnames_ta)
  deallocate(time_qa); deallocate(qa_in); deallocate(fnames_qa)
  deallocate(time_uw); deallocate(uw_in); deallocate(fnames_uw)
  deallocate(time_vw); deallocate(vw_in); deallocate(fnames_vw)
#ifdef WSPEED_IN
  deallocate(time_ws); deallocate(ws_in); deallocate(fnames_ws)
#endif
  deallocate(fnames_clm_temp);  deallocate(fnames_clm_salt)
  deallocate(fnames_clm_u);  deallocate(fnames_clm_v)
  deallocate(fnames_adv_temp);  deallocate(fnames_adv_salt)
  deallocate(fnames_adv_u);  deallocate(fnames_adv_v)
  deallocate(penet)
  deallocate(temp_adv_1d) ; deallocate(salt_adv_1d)
  deallocate(u_adv_1d);  deallocate(v_adv_1d)
end program main
