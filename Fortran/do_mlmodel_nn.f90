program do_mlmodel_nn
  ! Main program of solving Nakanishi-Niino ML model
  use ml_param
  use ml_utils
  use nnf_sub
  use solve_diag
  implicit none
  real(idx),parameter :: pi=4.0_idx*atan(1.0_idx),omega=7.29e-5_idx
  integer :: nz_rho,nz_q
  real(idx) :: lat,coriolis
  real(idx),allocatable :: z_rho(:),z_q(:)
  real(idx),allocatable :: lev_rho(:),lev_q(:)
  real(idx),allocatable :: temp_1d(:),salt_1d(:)
  real(idx),allocatable :: u_1d(:),v_1d(:)
  real(idx),allocatable :: temp_next_1d(:),salt_next_1d(:)
  real(idx),allocatable :: u_next_1d(:),v_next_1d(:)
  real(idx),allocatable :: temp_past_1d(:),salt_past_1d(:)
  real(idx),allocatable :: u_past_1d(:),v_past_1d(:)
  real(idx),allocatable :: temp_correct_1d(:),salt_correct_1d(:)
  real(idx),allocatable :: u_correct_1d(:),v_correct_1d(:)
  real(idx),allocatable :: kt_1d(:),ks_1d(:),km_1d(:)
  real(idx),allocatable :: qq_1d(:),qq_next_1d(:),qq_past_1d(:)
  real(idx),allocatable :: l_1d(:),l_next_1d(:)
  real(idx),allocatable :: kq_1d(:)
  real(idx),allocatable :: shear_1d(:),bvf_1d(:)
  real(idx),allocatable :: penet_1d(:)
  real(idx) :: hflx_solar,hflx_nosolar
  real(idx) :: sflx
  real(idx) :: uflx,vflx
  real(idx) :: bottom_depth,dz
  integer :: iz,itime,ntime,ntime_output,nstep_out
  real(idx) :: dt,dt_output
  type(datetime) :: dt_start,dt_end,dt_now
  real(idx) :: jd_start,jd_end,jd_now
  character(len=maxlen) :: fname_out
  ! Input parameter
  dt=600.0
  dt=240.0
  dt_start%year=1900;dt_start%month=1;dt_start%day=1
  dt_start%hour=0;dt_start%minute=0;dt_start%second=0
  dt_end%year=1900;dt_end%month=1;dt_end%day=21
  dt_output=60.0_idx*60.0_idx*24.0_idx
  dt_end%hour=0;dt_end%minute=0;dt_end%second=0

  bottom_depth=1000.0;dz=5
  lat=45.0
  fname_out="out_case_nn.txt"
  nu=1.0e-6     ! background momentum diffusion
  nu_t=1.0e-7   ! background temperature diffusion
  nu_s=1.0e-7   ! background salinity diffusion
  beta1=0.6     ! Shortwave penetration parameter (Jerov water type) [m]
  beta2=20.0    ! Shortwave penetration parameter [m]
  r_long=0.62   ! Shortwave penetration parameter [non dimensional]

 
  jd_start=dt_to_jd(dt_start)
  jd_end=dt_to_jd(dt_end)
  ntime=int((jd_end-jd_start)*day_to_sec/dt)
  ntime_output=int(dt_output/dt)

  coriolis=2.0*omega*sin(lat*pi/180.0)
  nz_rho=int(bottom_depth/dz)
  nz_q=nz_rho+1
  allocate(z_rho(1:nz_rho));  allocate(lev_rho(1:nz_rho))
  allocate(z_q(1:nz_q));  allocate(lev_q(1:nz_q))
  do iz = 1,nz_q
     z_q(iz)=-1.0_idx*(iz-1)*dz
     lev_q(iz)=-1.0_idx*z_q(iz)
  end do
  do iz=1,nz_rho
     z_rho(iz)=0.5_idx*(z_q(iz)+z_q(iz+1))
     lev_rho(iz)=0.5_idx*(lev_q(iz)+lev_q(iz+1))
  end do

  allocate(temp_1d(1:nz_rho));allocate(salt_1d(1:nz_rho))
  allocate(u_1d(nz_rho));allocate(v_1d(nz_rho))
  allocate(temp_next_1d(nz_rho));allocate(salt_next_1d(nz_rho))
  allocate(u_next_1d(nz_rho));allocate(v_next_1d(nz_rho))
  allocate(temp_past_1d(nz_rho));allocate(salt_past_1d(nz_rho))
  allocate(u_past_1d(nz_rho));allocate(v_past_1d(nz_rho))
  allocate(temp_correct_1d(nz_rho));allocate(salt_correct_1d(nz_rho))
  allocate(u_correct_1d(nz_rho));allocate(v_correct_1d(nz_rho))
  temp_correct_1d=0.0_idx;salt_correct_1d=0.0_idx
  u_correct_1d=0.0_idx;v_correct_1d=0.0_idx
  allocate(bvf_1d(1:(nz_rho-1))); allocate(shear_1d(1:(nz_rho-1)))
  allocate(km_1d(1:(nz_rho-1))); allocate(kt_1d(1:(nz_rho-1))) ; allocate(ks_1d(1:(nz_rho-1)))
  bvf_1d=0.0_idx; shear_1d=0.0_idx ; km_1d=0.0_idx; kt_1d=0.0_idx; ks_1d=0.0_idx
  allocate(penet_1d(nz_q))
  allocate(qq_1d(1:(nz_q)));allocate(qq_next_1d(1:(nz_q)))
  allocate(qq_past_1d(1:(nz_q)))
  allocate(l_1d(1:(nz_q)));allocate(l_next_1d(1:(nz_q)))
  allocate(kq_1d(1:(nz_rho)))
  kq_1d=0.0_idx
  open(10,file=fname_out,status="replace",action="write")

  !========================!
  ! Set initial condition  !
  !========================!
  do iz =1,nz_rho
     if (lev_rho(iz) .le. 20.0_idx) then
        temp_1d(iz)=15.0_idx
        salt_1d(iz)=35.0_idx
     else if (lev_rho(iz) .le. 50.0_idx) then
        temp_1d(iz)=15.0_idx+(12.0_idx-15.0_idx)*(lev_rho(iz)-20.0_idx)/(50.0_idx-20.0_idx)
     else
        temp_1d(iz)=12.0_idx
     end if
     salt_1d(iz)=35.0_idx
     u_1d(iz)=0.0_idx
     v_1d(iz)=0.0_idx
  end do
  call initialize_turb(nz_q,z_q,qq_1d,l_1d)
  temp_past_1d=temp_1d;salt_past_1d=salt_1d;u_past_1d=u_1d;v_past_1d=v_1d
  qq_past_1d=qq_1d
  nstep_out=0
  do itime = 1,ntime
     if (mod(itime,ntime_output) == 0) then
        nstep_out=nstep_out+1
     end if
  end do
  write(10,*) nstep_out
  write(10,*) (lev_rho(iz),iz=1,nz_rho)
  do itime = 1,ntime
     jd_now=jd_start+dt*(itime)/day_to_sec
     ! Case 1 wind deepening
     ! No heat/salt flux, only zonal wind stress with amplitude of 0.2 [N/m^2]
     hflx_nosolar=0.0_idx;hflx_solar=0.0_idx
     uflx=0.2;vflx=0.0
     !   Case 2
     !hflx_nosolar=-30.0_idx;hflx_solar=0.0_idx
     !sflx=0.0_idx
     !uflx=0.0;vflx=0.0
     !  Case 3
     !hflx_nosolar=30.0_idx;hflx_solar=0.0_idx
     !sflx=0.0_idx
     !  Case 4
     hflx_nosolar=-30.0_idx;hflx_solar=0.0_idx
     sflx=0.0_idx
     uflx=0.2;vflx=0.0
     !  Case 5
     sflx=0.0_idx
     uflx=0.2;vflx=0.0
     sflx=-2.5e-5 ! (E-P)
     ! Calculate shear and stratification
     bvf_1d=cal_bvf(nz_rho,z_rho,temp_1d,salt_1d)
     shear_1d=cal_shear(nz_rho,z_rho,u_1d,v_1d)
     call heat_absorb(nz_q,z_q,hflx_solar,beta1,beta2,r_long,penet_1d)
     call mynnf25_vdiff(nz_rho,bvf_1d,shear_1d,qq_1d,l_1d,&
                   & km_1d,kt_1d,ks_1d)
     do iz=1,nz_rho-1
       km_1d(iz) = km_1d(iz)  + nu ; kt_1d(iz) = kt_1d(iz)  + nu_t
       ks_1d(iz) = ks_1d(iz)  + nu_s
    end do
     call cal_kq(nz_rho,km_1d,kq_1d)
     
     call cal_next_theta(2*dt,nz_rho,z_rho,z_q,temp_past_1d,kt_1d,penet_1d,&
          & hflx_nosolar,temp_correct_1d,temp_next_1d)
     call cal_next_salt(2*dt,nz_rho,z_rho,z_q,salt_past_1d,&
          & ks_1d,sflx,salt_correct_1d,salt_next_1d)
     call cal_next_uv(2*dt,nz_rho,z_rho,z_q,u_past_1d,v_past_1d,&
          & km_1d,coriolis,uflx,vflx,u_correct_1d,v_correct_1d,&
          & u_next_1d,v_next_1d)
     call cal_next_qq(2*dt,nz_rho,z_rho,z_q,qq_past_1d,bvf_1d,shear_1d, &
          & l_1d,kq_1d,km_1d,kt_1d,uflx,vflx,B_1,qq_next_1d)
     call diag_l_mynnf(nz_rho,z_q,bvf_1d,qq_1d,l_next_1d,&
          & uflx,vflx,hflx_solar,hflx_nosolar,sflx)
     !================!
     ! Update values  !
     !================!
     do iz =1,nz_rho
        temp_past_1d(iz)=temp_1d(iz) ;  salt_past_1d(iz)=salt_1d(iz)
        temp_1d(iz)=temp_next_1d(iz) ;  salt_1d(iz)=salt_next_1d(iz)
        u_past_1d(iz)=u_1d(iz) ;  v_past_1d(iz)=v_1d(iz)
        u_1d(iz)=u_next_1d(iz) ;  v_1d(iz)=v_next_1d(iz)
     end do
     do iz=1,nz_q
        qq_past_1d(iz) = qq_1d(iz)
        qq_1d(iz) = qq_next_1d(iz)
        l_1d(iz) = l_next_1d(iz)
     end do
     if (mod(itime,ntime_output) == 0) then
        write(*,*) dz*dz/(4*maxval(km_1d)),dt
        dt_now=jd_to_dt(jd_now)
        write(*,*) dt_now,temp_1d(1)
        write(10,'(6i6)') dt_now%year,dt_now%month,dt_now%day,&
             & dt_now%hour,dt_now%minute,dt_now%second
        write(10,*) (temp_1d(iz),iz=1,nz_rho)
        write(10,*) (salt_1d(iz),iz=1,nz_rho)
        write(10,*) (u_1d(iz),iz=1,nz_rho)
        write(10,*) (v_1d(iz),iz=1,nz_rho)
     end if
     !   call
  end do
  close(10)
  deallocate(z_rho);deallocate(z_q)
  deallocate(lev_rho);deallocate(lev_q)
  deallocate(temp_1d);deallocate(salt_1d);deallocate(u_1d);deallocate(v_1d)
  deallocate(temp_next_1d);deallocate(salt_next_1d)
  deallocate(u_next_1d);deallocate(v_next_1d)
  deallocate(temp_past_1d);deallocate(salt_past_1d)
  deallocate(u_past_1d);deallocate(v_past_1d)
  deallocate(temp_correct_1d);deallocate(salt_correct_1d)
  deallocate(u_correct_1d);deallocate(v_correct_1d)
  deallocate(kt_1d);deallocate(ks_1d);deallocate(km_1d)
  deallocate(bvf_1d);deallocate(shear_1d)
  deallocate(penet_1d)
  deallocate(qq_1d);deallocate(qq_next_1d);deallocate(qq_past_1d)
  deallocate(l_1d);deallocate(l_next_1d)
  deallocate(kq_1d)
end program do_mlmodel_nn
