program lin_main
  ! Main source program of mixed layer model
  use param
  use arrays_sub
  use input_files
  use mod_adv
  use mod_atm
  use mod_hist
  use mod_diags
  use mod_lin
!  use mod_rst
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
  character(maxlen) :: fname_grid,fname_init_rho,fname_init_q
  integer :: nlon,nlat,nz,ntime
  ! Time step parameters
  real(idx) :: total_time,elapsed_time,tmp
  real(idx),parameter :: eps=5.0e-3
  real(idx),allocatable :: temp_ori(:,:,:),salt_ori(:,:,:),u_ori(:,:,:),v_ori(:,:,:)
  real(idx),allocatable :: qq_ori(:,:,:),l_ori(:,:,:)
  real(idx),allocatable :: temp_rate_ori(:,:,:),salt_rate_ori(:,:,:)
  real(idx),allocatable :: u_rate_ori(:,:,:),v_rate_ori(:,:,:)
!  real(idx),allocatable :: matrix(:,:,:,:)
  !===========================================
  ! Output setting
  integer :: ix,iy,iz,it,imat
  !$ double precision st, en
  !===========================================
  namelist/input_init/fname_grid,fname_init_rho,fname_init_q
  !$ st = omp_get_wtime()
  !===========================================!
  ! Read namelists                            !
  !===========================================!
  open(unit=nmlf_master,file=trim(namelist_master))
  read(unit=nmlf_master,nml=master)
  !===========================================!
  ! Load parameters                           !
  !===========================================!
  open(unit=nmlf_set,file=trim(namelist_set))  
  read(unit=nmlf_set,nml=parameters)
  close(nmlf_set)
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
  read(unit=nmlf_io,nml=input_init)
  call read_grid_file(fname_grid,nlon,nlat,nz,lon_grd,lat_grd,cor,z_rho,z_q)
  call allocate_ocn_3d_arrays(nlon,nlat,nz,temp,salt,u,v,qq,l)
  call allocate_ocn_3d_arrays(nlon,nlat,nz,temp_ori,salt_ori,u_ori,v_ori,qq_ori,l_ori)
  call allocate_atm_2d_arrays(nlon,nlat,sw,lw,sh,lh,ev,pr,ust,vst)
  !===========================================!
  ! Read initial contidion                    !
  !===========================================!
  call read_init_rho_file(fname_init_rho,nlon,nlat,nz,temp,salt,u,v)
  call read_init_rho_file(fname_init_rho,nlon,nlat,nz,temp_ori,salt_ori,u_ori,v_ori)
  ! Initialize turbulent field-----------------------------------------------
#if defined(NNF) | defined(KC) | defined(bd)
  qq = 0.0_idx ; l=0.0_idx
  call read_init_q_file(fname_init_q,nlon,nlat,nz,qq,l)
  call read_init_q_file(fname_init_q,nlon,nlat,nz,qq_ori,l_ori)
#endif
  !===========================================!
  ! Read atmospheric forcing file             !
  !===========================================!
  read(unit=nmlf_io,nml=input_atm)
  ! Precipitation
  call read_atmos_file(fname_pr,varname_pr,ntime_pr,time_pr,pr_in,start_yymmdd_int,start_hhmmss_int)
  ! Radiation
  call read_atmos_file(fname_sw,varname_sw,ntime_sw,time_sw,sw_in,start_yymmdd_int,start_hhmmss_int)
  call read_atmos_file(fname_lw,varname_lw,ntime_lw,time_lw,lw_in,start_yymmdd_int,start_hhmmss_int)
  call read_atmos_file(fname_ta,varname_ta,ntime_ta,time_ta,ta_in,start_yymmdd_int,start_hhmmss_int)
  call read_atmos_file(fname_qa,varname_qa,ntime_qa,time_qa,qa_in,start_yymmdd_int,start_hhmmss_int)
  call read_atmos_file(fname_uw,varname_uw,ntime_uw,time_uw,uw_in,start_yymmdd_int,start_hhmmss_int)
  call read_atmos_file(fname_vw,varname_vw,ntime_vw,time_vw,vw_in,start_yymmdd_int,start_hhmmss_int)
#if defined(WSPEED_IN)
  call read_atmos_file(fname_ws,varname_ws,ntime_ws,time_ws,ws_in,start_yymmdd_int,start_hhmmss_int)  
#endif
  write(*,*) "*******************************************************"
  write(*,*) " Finish reading atmospheric file"
  write(*,*) "*******************************************************"
  !===========================================!
  ! Read advection file                       !
  !===========================================!
  read(unit=nmlf_io,nml=input_adv)
  if (use_adv_temp .eqv. .TRUE.) then
     call read_adv_file(fname_adv_temp,varname_adv_temp,ntime_adv_temp,time_adv_temp,&
          & temp_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  if (use_adv_salt .eqv. .TRUE.) then
     call read_adv_file(fname_adv_salt,varname_adv_salt,ntime_adv_salt,time_adv_salt,&
          & salt_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  if (use_adv_u .eqv. .TRUE.) then
     call read_adv_file(fname_adv_u,varname_adv_u,ntime_adv_u,time_adv_u,&
          & u_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  if (use_adv_v .eqv. .TRUE.) then
     call read_adv_file(fname_adv_v,varname_adv_v,ntime_adv_v,time_adv_v,&
          & v_adv,start_yymmdd_int,start_hhmmss_int)
  end if
  !=================================================
  ! Allocate arrays necessary for calculation
  !==================================================
  call allocate_mix_3d_arrays(nlon,nlat,nz,bvf,shear,km,kt)
  call allocate_ocn_1d_arrays(nz,temp_next,salt_next,u_next,v_next,qq_next,l_next)
  call allocate_mix_1d_arrays(nz-1,bvf_1d,shear_1d,km_1d,kt_1d,ks_1d,kq_1d)
  allocate(penet(nz))
  allocate(temp_adv_1d(nz)) ; allocate(salt_adv_1d(nz))
  allocate(u_adv_1d(nz)) ;   allocate(v_adv_1d(nz))
  allocate(matrix(nlon,nlat,4*nz,4*nz))

  allocate(temp_rate_ori(1:nlon,1:nlat,1:nz))
  allocate(salt_rate_ori(1:nlon,1:nlat,1:nz))
  allocate(u_rate_ori(1:nlon,1:nlat,1:nz))
  allocate(v_rate_ori(1:nlon,1:nlat,1:nz))
  temp_next=0.0_idx ; salt_next =0.0_idx;  u_next = 0.0_idx  ; v_next = 0.0_idx
  qq_next = 0.0_idx ; l_next=0.0_idx
  penet=0.0_idx
  !======================!
  ! Prepare output files !
  !======================!
  !*********************************
  ! Diagnostic (temperature) file
  !*********************************
  read(unit=nmlf_io,nml=output_temp_diag_flag)
  call allocate_tdiag(nlon,nlat,nz,temp_rate,temp_vdiff,temp_penet,temp_relax)
  !*********************************
  ! Diagnostic (salinity) file
  !*********************************
  read(unit=nmlf_io,nml=output_salt_diag_flag)
  call allocate_sdiag(nlon,nlat,nz,salt_rate,salt_vdiff,salt_relax)
  !*********************************
  ! Diagnostic (uv) file
  !*********************************
  read(unit=nmlf_io,nml=output_uv_diag_flag)
  call allocate_uvdiag(nlon,nlat,nz,u_rate,u_cor,u_vdiff,u_relax,v_rate,v_cor,v_vdiff,v_relax)
  !*********************************
  ! Diagnostic (qq) file
  !*********************************
  read(unit=nmlf_io,nml=output_qq_diag_flag)
  call allocate_qqdiag(nlon,nlat,nz,qq_rate,qq_vdiff,qq_sp,qq_bp,qq_disp)
  close(nmlf_io)
  !==========================!
  !  Initialize output flag  !
  !==========================!
  !==========================!
  !  Start point of loop     !
  !==========================!
  ntime=1
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
  !$omp parallel do private(ld,ta,qa,uwind,vwind,ws,heat_solar,heat_no_solar,e_p,&
  !$omp & penet,bvf_1d,shear_1d,km_1d,kt_1d,ks_1d,kq_1d,&
  !$omp & temp_adv_1d,salt_adv_1d,u_adv_1d,v_adv_1d,&
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
           if (use_adv_temp .eqv. .TRUE.) then
              do iz = 1,nz
                 temp_adv_1d(iz)=set_data(ind1_adv_temp,ind2_adv_temp,wgt1_adv_temp,wgt2_adv_temp,&
                      & temp_adv(ix,iy,iz,:))
              end do
           else                
              temp_adv_1d = 0.0_idx
           end if
           ! Salinity equation
           if (use_adv_salt .eqv. .TRUE.) then
              do iz = 1,nz
                 salt_adv_1d(iz)=set_data(ind1_adv_salt,ind2_adv_salt,wgt1_adv_salt,wgt2_adv_salt,&
                      & salt_adv(ix,iy,iz,:))
              end do
           else
              salt_adv_1d = 0.0_idx                
           end if
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
           !              call cal_next_qq(dt,nz,z_rho,z_q,qq(ix,iy,:),bvf_1d(1:nz-1),shear_1d(1:nz-1) &
           !                   & ,l(ix,iy,:),kq_1d,km_1d,kt_1d,ust(ix,iy),vst(ix,iy),B_1,qq_next)
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
           call cal_next_theta_bud(dt,nz,z_rho,z_q,temp(ix,iy,1:nz),kt_1d,penet,&
                & heat_no_solar,temp_adv_1d,temp_next,&
                & temp_rate_ori(ix,iy,1:nz),temp_vdiff(ix,iy,1:nz),&
                & temp_penet(ix,iy,1:nz),temp_relax(ix,iy,1:nz))
           call cal_next_salt_bud(dt,nz,z_rho,z_q,salt(ix,iy,1:nz),&
                & ks_1d,ssflux,salt_adv_1d,salt_next,&
                & salt_rate_ori(ix,iy,1:nz),salt_vdiff(ix,iy,1:nz),salt_relax(ix,iy,1:nz))
           call cal_next_uv_bud(dt,nz,z_rho,z_q,u(ix,iy,:),v(ix,iy,:),&
                & km_1d,cor(iy),ust(ix,iy),vst(ix,iy),u_adv_1d,v_adv_1d,u_next,v_next,&
                & u_rate_ori(ix,iy,1:nz),u_cor(ix,iy,1:nz),u_vdiff(ix,iy,1:nz),u_relax(ix,iy,1:nz),&
                & v_rate_ori(ix,iy,1:nz),v_cor(ix,iy,1:nz),v_vdiff(ix,iy,1:nz),v_relax(ix,iy,1:nz))

           ! Temperature perturbation
           do imat = 1,nz
              temp=temp_ori ; salt=salt_ori; u=u_ori; v=v_ori
              qq=qq_ori ; l=l_ori
              temp_rate(ix,iy,1:nz)=0.0_idx; salt_rate(ix,iy,1:nz)=0.0_idx
              u_rate(ix,iy,1:nz)=0.0_idx; v_rate(ix,iy,1:nz)=0.0_idx
              temp(ix,iy,imat)=temp_ori(ix,iy,imat)+eps
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
              call cal_next_theta_bud(dt,nz,z_rho,z_q,temp(ix,iy,1:nz),kt_1d,penet,&
                   & heat_no_solar,temp_adv_1d,temp_next,&
                   & temp_rate(ix,iy,1:nz),temp_vdiff(ix,iy,1:nz),&
                   & temp_penet(ix,iy,1:nz),temp_relax(ix,iy,1:nz))
              call cal_next_salt_bud(dt,nz,z_rho,z_q,salt(ix,iy,1:nz),&
                   & ks_1d,ssflux,salt_adv_1d,salt_next,&
                   & salt_rate(ix,iy,1:nz),salt_vdiff(ix,iy,1:nz),salt_relax(ix,iy,1:nz))
              call cal_next_uv_bud(dt,nz,z_rho,z_q,u(ix,iy,:),v(ix,iy,:),&
                   & km_1d,cor(iy),ust(ix,iy),vst(ix,iy),u_adv_1d,v_adv_1d,u_next,v_next,&
                   & u_rate(ix,iy,1:nz),u_cor(ix,iy,1:nz),u_vdiff(ix,iy,1:nz),u_relax(ix,iy,1:nz),&
                   & v_rate(ix,iy,1:nz),v_cor(ix,iy,1:nz),v_vdiff(ix,iy,1:nz),v_relax(ix,iy,1:nz))
              matrix(ix,iy,1:nz,imat)=(temp_rate(ix,iy,1:nz)-temp_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,nz+1:2*nz,imat)=(salt_rate(ix,iy,1:nz)-salt_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,2*nz+1:3*nz,imat)=(u_rate(ix,iy,1:nz)-u_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,3*nz+1:4*nz,imat)=(v_rate(ix,iy,1:nz)-v_rate_ori(ix,iy,1:nz))/eps
           end do
           ! Salinity perturbation
          do imat = 1,nz
             temp=temp_ori ; salt=salt_ori; u=u_ori; v=v_ori
             qq=qq_ori ; l=l_ori
             temp_rate(ix,iy,1:nz)=0.0_idx; salt_rate(ix,iy,1:nz)=0.0_idx
             u_rate(ix,iy,1:nz)=0.0_idx; v_rate(ix,iy,1:nz)=0.0_idx
             salt(ix,iy,imat)=salt_ori(ix,iy,imat)+eps
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
              call cal_next_theta_bud(dt,nz,z_rho,z_q,temp(ix,iy,1:nz),kt_1d,penet,&
                   & heat_no_solar,temp_adv_1d,temp_next,&
                   & temp_rate(ix,iy,1:nz),temp_vdiff(ix,iy,1:nz),&
                   & temp_penet(ix,iy,1:nz),temp_relax(ix,iy,1:nz))
              call cal_next_salt_bud(dt,nz,z_rho,z_q,salt(ix,iy,1:nz),&
                   & ks_1d,ssflux,salt_adv_1d,salt_next,&
                   & salt_rate(ix,iy,1:nz),salt_vdiff(ix,iy,1:nz),salt_relax(ix,iy,1:nz))
              call cal_next_uv_bud(dt,nz,z_rho,z_q,u(ix,iy,:),v(ix,iy,:),&
                   & km_1d,cor(iy),ust(ix,iy),vst(ix,iy),u_adv_1d,v_adv_1d,u_next,v_next,&
                   & u_rate(ix,iy,1:nz),u_cor(ix,iy,1:nz),u_vdiff(ix,iy,1:nz),u_relax(ix,iy,1:nz),&
                   & v_rate(ix,iy,1:nz),v_cor(ix,iy,1:nz),v_vdiff(ix,iy,1:nz),v_relax(ix,iy,1:nz))
              matrix(ix,iy,1:nz,imat+nz*1)=(temp_rate(ix,iy,1:nz)-temp_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,nz+1:2*nz,imat+nz*1)=(salt_rate(ix,iy,1:nz)-salt_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,2*nz+1:3*nz,imat+nz*1)=(u_rate(ix,iy,1:nz)-u_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,3*nz+1:4*nz,imat+nz*1)=(v_rate(ix,iy,1:nz)-v_rate_ori(ix,iy,1:nz))/eps
           end do

           ! U perturbation
           do imat = 1,nz
              temp=temp_ori ; salt=salt_ori; u=u_ori; v=v_ori
              qq=qq_ori ; l=l_ori
              temp_rate(ix,iy,1:nz)=0.0_idx; salt_rate(ix,iy,1:nz)=0.0_idx
              u_rate(ix,iy,1:nz)=0.0_idx; v_rate(ix,iy,1:nz)=0.0_idx
              u(ix,iy,imat)=u_ori(ix,iy,imat)+eps
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
              call cal_next_theta_bud(dt,nz,z_rho,z_q,temp(ix,iy,1:nz),kt_1d,penet,&
                   & heat_no_solar,temp_adv_1d,temp_next,&
                   & temp_rate(ix,iy,1:nz),temp_vdiff(ix,iy,1:nz),&
                   & temp_penet(ix,iy,1:nz),temp_relax(ix,iy,1:nz))
              call cal_next_salt_bud(dt,nz,z_rho,z_q,salt(ix,iy,1:nz),&
                   & ks_1d,ssflux,salt_adv_1d,salt_next,&
                   & salt_rate(ix,iy,1:nz),salt_vdiff(ix,iy,1:nz),salt_relax(ix,iy,1:nz))
              call cal_next_uv_bud(dt,nz,z_rho,z_q,u(ix,iy,:),v(ix,iy,:),&
                   & km_1d,cor(iy),ust(ix,iy),vst(ix,iy),u_adv_1d,v_adv_1d,u_next,v_next,&
                   & u_rate(ix,iy,1:nz),u_cor(ix,iy,1:nz),u_vdiff(ix,iy,1:nz),u_relax(ix,iy,1:nz),&
                   & v_rate(ix,iy,1:nz),v_cor(ix,iy,1:nz),v_vdiff(ix,iy,1:nz),v_relax(ix,iy,1:nz))
              matrix(ix,iy,1:nz,imat+nz*2)=(temp_rate(ix,iy,1:nz)-temp_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,nz+1:2*nz,imat+nz*2)=(salt_rate(ix,iy,1:nz)-salt_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,2*nz+1:3*nz,imat+nz*2)=(u_rate(ix,iy,1:nz)-u_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,3*nz+1:4*nz,imat+nz*2)=(v_rate(ix,iy,1:nz)-v_rate_ori(ix,iy,1:nz))/eps
           end do

           ! V perturbation
           do imat = 1,nz
              temp=temp_ori ; salt=salt_ori; u=u_ori; v=v_ori
              qq=qq_ori ; l=l_ori
              temp_rate(ix,iy,1:nz)=0.0_idx; salt_rate(ix,iy,1:nz)=0.0_idx
              u_rate(ix,iy,1:nz)=0.0_idx; v_rate(ix,iy,1:nz)=0.0_idx
              v(ix,iy,imat)=v_ori(ix,iy,imat)+eps
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
              !              call cal_next_qq(dt,nz,z_rho,z_q,qq(ix,iy,:),bvf_1d(1:nz-1),shear_1d(1:nz-1) &
              !                   & ,l(ix,iy,:),kq_1d,km_1d,kt_1d,ust(ix,iy),vst(ix,iy),B_1,qq_next)
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
              call cal_next_theta_bud(dt,nz,z_rho,z_q,temp(ix,iy,1:nz),kt_1d,penet,&
                   & heat_no_solar,temp_adv_1d,temp_next,&
                   & temp_rate(ix,iy,1:nz),temp_vdiff(ix,iy,1:nz),&
                   & temp_penet(ix,iy,1:nz),temp_relax(ix,iy,1:nz))
              call cal_next_salt_bud(dt,nz,z_rho,z_q,salt(ix,iy,1:nz),&
                   & ks_1d,ssflux,salt_adv_1d,salt_next,&
                   & salt_rate(ix,iy,1:nz),salt_vdiff(ix,iy,1:nz),salt_relax(ix,iy,1:nz))
              call cal_next_uv_bud(dt,nz,z_rho,z_q,u(ix,iy,:),v(ix,iy,:),&
                   & km_1d,cor(iy),ust(ix,iy),vst(ix,iy),u_adv_1d,v_adv_1d,u_next,v_next,&
                   & u_rate(ix,iy,1:nz),u_cor(ix,iy,1:nz),u_vdiff(ix,iy,1:nz),u_relax(ix,iy,1:nz),&
                   & v_rate(ix,iy,1:nz),v_cor(ix,iy,1:nz),v_vdiff(ix,iy,1:nz),v_relax(ix,iy,1:nz))
              matrix(ix,iy,1:nz,imat+nz*3)=(temp_rate(ix,iy,1:nz)-temp_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,nz+1:2*nz,imat+nz*3)=(salt_rate(ix,iy,1:nz)-salt_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,2*nz+1:3*nz,imat+nz*3)=(u_rate(ix,iy,1:nz)-u_rate_ori(ix,iy,1:nz))/eps
              matrix(ix,iy,3*nz+1:4*nz,imat+nz*3)=(v_rate(ix,iy,1:nz)-v_rate_ori(ix,iy,1:nz))/eps
           end do
           
!           do imat=1,nz
!              write(*,*) imat
!              write(*,*) matrix(nz+imat,1:4*nz)
 !             write(*,*) 
 !          end do
        else
           !    temp(ix,iy,:) = missing_value ; salt(ix,iy,:) = missing_value
           !    u(ix,iy,:) = missing_value ; v(ix,iy,:) = missing_value
           !    bvf(ix,iy,1:nz-1)=missing_value
           !    shear(ix,iy,1:nz-1)=missing_value
           !    km(ix,iy,1:nz-1)=missing_value
           !    kt(ix,iy,1:nz-1)=missing_value
        end if
        write(*,*) matrix(ix,iy,30,35)
     end do
  end do
  !==================!
  ! output procedure !
  !==================!
  open(unit=nmlf_io,file=trim(namelist_io))
  read(unit=nmlf_io,nml=output_lin_flag)
  call prepare_output_lin(fname_out_lin_rho,nlon,nlat,nz,lon_grd,lat_grd,z_rho,missing_value)
  call output_lin_rho(fname_out_lin_rho,nlon,nlat,nz,matrix)
  close(nmlf_io)
  !      ! Deallocation
  deallocate(lon_grd); deallocate(lat_grd); deallocate(z_rho)
  deallocate(z_q)
  call deallocate_ocn_3d_arrays(temp,salt,u,v,qq,l)
  call deallocate_atm_2d_arrays(sw,lw,sh,lh,ev,pr,ust,vst)
  call deallocate_mix_3d_arrays(bvf,shear,km,kt)
  call deallocate_ocn_1d_arrays(temp_next,salt_next,u_next,v_next,qq_next,l_next)
  call deallocate_mix_1d_arrays(bvf_1d,shear_1d,km_1d,kt_1d,ks_1d,kq_1d)
  call deallocate_tdiag(temp_rate,temp_vdiff,temp_penet,temp_relax)
  !      ! Atmospheric data
  deallocate(time_pr); deallocate(pr_in)
  deallocate(time_sw); deallocate(sw_in)
  deallocate(time_lw); deallocate(lw_in)
  deallocate(time_ta); deallocate(ta_in)
  deallocate(time_qa); deallocate(qa_in)
  deallocate(time_uw); deallocate(uw_in)
  deallocate(time_vw); deallocate(vw_in)
  ! #ifdef WSPEED_IN
  !      deallocate(time_ws); deallocate(ws_in)
  ! #endif
  !      deallocate(penet)
  !      deallocate(temp_adv_1d) ; deallocate(salt_adv_1d)
  !      deallocate(u_adv_1d);  deallocate(v_adv_1d)
end program lin_main
