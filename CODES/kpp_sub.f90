module kpp_sub
  use param
  use eos_sub
  implicit none
  private
  !------------------------------------------
  real(idx),parameter :: epsilon = 0.1_idx,betaT=-0.2_idx
  ! Richardson mixing
  real(idx),parameter :: Ri0 = 0.7_idx
  real(idx),parameter :: anu0=50.0e-4_idx,anum=1.0e-4_idx,anut=1.0e-5_idx
  ! Surface
  real(idx),parameter :: zetas=-1.0_idx,zetam=-0.2_idx
  real(idx),parameter :: Cv=1.5_idx,Ric=0.30_idx
  real(idx),parameter :: tiny=1.0e-20_idx
  
  real(idx),parameter :: cs = 24.0_idx * (1.0_idx-16.0_idx*zetas)**(0.5_idx)
  real(idx),parameter :: cm = 12.0_idx * (1.0_idx-16.0_idx*zetam)**(-0.25_idx)
  real(idx),parameter :: as = cs * zetas + (1.0_idx -16.0_idx *zetas) **  (1.5_idx)
  real(idx),parameter :: am = cm * zetam + (1.0_idx -16.0_idx *zetam) **  (0.75_idx)
  real(idx),parameter :: Vtc = Cv*sqrt(-betaT)/(sqrt(cs*epsilon)*Ric*kappa*kappa)
  ! public setting--------------------------------------------------------
  public :: cal_hbl
  public :: cal_mo_inv,kpp_vdiff
contains
  !============================================================================
  function heat_to_bf(heat) result(bf)
    implicit none
    real(idx),intent(in) :: heat
    real(idx) :: bf
    real(idx),parameter :: sst=20.0_idx,sss=35.0_idx
    real(idx) :: alpha
    alpha = cal_alpha(sss,sst) ! <0
    bf = -1.0_idx * g * alpha * heat / (rho * cpw)
  end function heat_to_bf
  ! Evap - precipitation to Buoyancy flux
  function ep_to_bf(ep) result(bf)
    implicit none
    real(idx),intent(in) :: ep
    real(idx) :: bf
    real(idx) :: beta
    real(idx),parameter :: sst=20.0_idx,sss=35.0_idx
    beta = cal_beta(sss,sst) ! <0
    bf = -1.0_idx * g * beta * ep * S0
  end function ep_to_bf
  !===========================================================================
  function cal_mo_inv(sst,sss,tau_x,tau_y,sw,nsw,e_p) result(mo_inv)
    real(idx),intent(in) :: sst,sss,tau_x,tau_y,sw,nsw,e_p
    real(idx) :: mo_inv
    real(idx) :: alpha,beta,bf,tau,u_star2,u_star3
    alpha=cal_alpha(sss,sst) ! < 0
    beta=cal_beta(sss,sst)   ! > 0
    !bf=-wb0=-g*alpha*wto+beta*ws0
    bf = -1.0_idx * g * (alpha * (sw+nsw) / (rho*cpw)+beta*S0*e_p)
    ! tau_x = rho * u_star^2 * (direction)
    ! tau_y = rho * v_star^2
    ! u_star3 = (sqrt(u_star^2+v_star^2))^3
    u_star2=sqrt(tau_x**2+tau_y**2)/rho
    u_star3 = sqrt(u_star2)**3
    !    mo =u_star3 / (dkappa * bf)
    u_star3=max(u_star3,1.0e-10)
    mo_inv = kappa * bf / u_star3
    !bf <0 -> mo < 0 ! unstable
    !bf >0 -> mo > 0 ! stable
  end function cal_mo_inv
  !========================================================================
  function phi_s(zeta) result(phi)
    implicit none
    real(idx),intent(in) :: zeta
    real(idx) :: phi
    if (zeta >= 0.0_idx) then
       phi = 1.0_idx + 5.0_idx * zeta
    else if (zeta > zetas) then
       phi = (1.0_idx - 16.0_idx * zeta)**(-1.0_idx / 2.0_idx)
    else
       phi = (as-cs*zeta)**(-1.0_idx / 3.0_idx)
    end if
  end function phi_s
  function phi_m(zeta) result(phi)
    implicit none
    real(idx),intent(in) :: zeta
    real(idx) :: phi
    if (zeta >= 0.0_idx) then
       phi = 1.0_idx + 5.0_idx * zeta
    else if (zeta > zetam) then
       phi = (1.0_idx - 16.0_idx * zeta)**(-1.0_idx / 4.0_idx)
    else
       phi = (am-cm*zeta)**(-1.0_idx / 3.0_idx)
    end if
  end function phi_m
  !========================================================================
  function cal_hbl(N,z_rho,z_q,T,S,U,V,taux,tauy,heat_solar,heat_no_solar,&
       & e_p,beta1,beta2,r_long,f) result(hbl)
    implicit none
    integer,intent(in) :: N
    real(idx),intent(in) :: z_rho(N), z_q(N)
    real(idx),intent(in) :: T(N),S(N),U(N),V(N)
    real(idx),intent(in) :: taux,tauy,heat_solar,heat_no_solar
    real(idx),intent(in) :: e_p,beta1,beta2,r_long,f
    real(idx) :: dens(N),Rib(N)
    integer :: iz,iter
    real(idx) :: dens_surf,mld,sl_depth
    real(idx) :: dens_ref,U_ref,V_ref
    real(idx) :: hbl,hbln,hmonob,hekman
    real(idx) :: ustar,Bo_net,Bo_solar,Bo_nonsolar
    real(idx) :: zeta,phis,w_s
    real(idx) :: delU,delV,dV2
    real(idx) :: bvtop,bvf2_1,dVt2,L_inv
    integer :: sl_ind
    integer,parameter :: maxiter=10
    do iz = 1,N
       dens(iz) = cal_dens(S(iz),T(iz))
    end do
    dens_surf=dens(N)+0.1_idx
    mld=abs(maxval(z_rho,MASK=(dens>=dens_surf)))
    hbl=mld
    ! Atmospheric forcing
    ustar=sqrt(sqrt(taux**2+tauy**2)/rho)
    Bo_nonsolar = heat_to_bf(heat_no_solar) + ep_to_bf(e_p)
              
    do iter=1,maxiter
       sl_depth=epsilon * hbl ; sl_depth=max(abs(z_rho(N)),sl_depth)
       sl_ind= sum(maxloc(z_rho,MASK=(abs(z_rho)>=sl_depth)))
       dens_ref = sum(dens(sl_ind:N)) / (N-sl_ind+1)
       u_ref= sum(u(sl_ind:N)) / (N-sl_ind+1); v_ref= sum(v(sl_ind:N)) / (N-sl_ind+1)

       do iz = 2,N-2
          Bo_solar = heat_to_bf(heat_solar - cal_penet(z_rho(iz),heat_solar,beta1,beta2,r_long))
          Bo_net = Bo_nonsolar+Bo_solar
          L_inv = kappa * Bo_net / (ustar**3+tiny)
          ! Velocity scales
          if (L_inv .ge. 0.0_idx) then
             zeta = abs(z_rho(iz)) * L_inv
          else
             zeta=min(abs(z_rho(iz)),sl_depth) * L_inv      
          end if
          phis = phi_s(zeta) ;   w_s = (kappa * ustar) / phis
          ! Delta
          delU = u_ref-u(iz) ;  delV = v_ref - v(iz)
          dV2 = delU*delU + delV*delV
          bvtop = -1.0_idx * (g / rho) * (dens_ref-dens(iz)) * abs(z_rho(iz))
          bvf2_1 = -1.0_idx * (g / rho) * (dens(iz)-dens(iz-1)) / (z_rho(iz)-z_rho(iz-1))
          dVt2 = Vtc*abs(z_rho(iz))*w_s*sqrt(max(bvf2_1,0.0_idx))
          Rib(iz) = bvtop / (dV2+dVt2+tiny)
       end do
       Rib(N)=Rib(N-1)
       Rib(1)=Rib(2)
       hbln=abs(maxval(z_rho,MASK=(Rib > Ric)))
       if (hbln .ne. hbl) then
          hbl = hbln
       else
          exit
       end if
    end do
    Bo_solar = heat_to_bf(heat_solar - cal_penet(-1.0*hbl,heat_solar,beta1,beta2,r_long))
    Bo_net = Bo_nonsolar+Bo_solar
    if (Bo_net .gt. 0) then
       hekman = 0.7_idx * ustar / (f+tiny)
       hmonob = ustar**3 / max(kappa * Bo_net,tiny)
       hbl = min(hekman,hmonob,hbl)
    end if
  end function cal_hbl
  !============================================================================
  subroutine kpp_vdiff(n,z_rho,z_q,temp,salt,u,v,&
       & bvf,shear2,KM,KT,KS,&
       & taux,tauy,heat_solar,heat_no_solar,e_p,beta1,beta2,r_long,f)
    implicit none
    integer,intent(in) :: n
    real(idx),intent(in) :: z_rho(N),z_q(N)
    real(idx) :: temp(n),salt(n),u(n),v(n)
    real(idx),intent(in) :: bvf(N-1),shear2(N-1)
    real(idx),intent(in) :: taux,tauy,heat_solar,heat_no_solar
    real(idx),intent(in) :: e_p,beta1,beta2,r_long,f
    real(idx),intent(inout) :: KM(N-1),KT(N-1),KS(N-1)
    real(idx) :: Rig(N-1),nu_ri(N-1)
    real(idx) :: ustar,Bo_solar,Bo_nonsolar,Bo_net
    real(idx) :: hbl
    real(idx) :: zlmd,sl_depth
    real(idx) :: L_inv,zetapar,phim,w_m,phis,w_s,f1
    real(idx) :: cff,cff_up,cff_dn,zeta,zetapar0
    real(idx) :: KMp,KTp,KMh,KTh,Gm1,Gt1,dGm1ds,dGt1ds
    real(idx) :: a2_m,a2_s,a3_m,a3_s
    real(idx) :: sigma,Gm,Gs
    integer :: i,k
    ! ================================================
    ! Interior mixing
    ! ================================================
    ! Richardson number dependent mixing
    Rig = bvf / (shear2 + tiny)
    do i = 1,N-1
       if (Rig(i) .lt. 0.0_idx) then
          nu_ri(i) = anu0 * 1.0_idx
       else if (Rig(i) .lt. Ri0) then
          nu_ri(i) = anu0 * (1.0_idx - (Rig(i)/Ri0)**2)**3
       else
          nu_ri(i) = 0.0_idx
       end if
    end do
    KM = nu_ri + anum ; KT = nu_ri + anut ; KS = nu_ri + anut

    ! Atmospheric forcing setting
    hbl=cal_hbl(n,z_rho,z_q,temp,salt,u,v,&
         & taux,tauy,heat_solar,heat_no_solar,e_p,beta1,beta2,r_long,f)
    ustar=sqrt(sqrt(taux**2+tauy**2)/rho)
    Bo_nonsolar = heat_to_bf(heat_no_solar) + ep_to_bf(e_p)
    Bo_solar = heat_to_bf(heat_solar - cal_penet(hbl,heat_solar,beta1,beta2,r_long))
    Bo_net = Bo_solar + Bo_nonsolar

    ! Interpolate the value at the bottom of hbl================
    !  Compute turbulent velocity scales (w_m,w_s) at hbl
    L_inv = kappa * Bo_net / (ustar**3+tiny)
    if (Bo_net >= 0) then
       zetapar  = 1.0_idx * hbl * L_inv
    else
       zetapar  = 1.0_idx * epsilon * hbl * L_inv       
    end if
    phim = phi_m(zetapar) ; phis = phi_s(zetapar)
    w_m = kappa * ustar / phim ; w_s = kappa * ustar / phis 
    ! Shape functions and derivatives at hbl( See LMD equation (18))
    if (Bo_net >= 0) then
       f1 = 5.0_idx * Bo_net * kappa / ((ustar**4)+tiny)
    else
       f1 = 0.0_idx
    end if

    ! Interpolate Km, Ks among grid points surrounding hbl
    ! compute G and derivative
    ! Situation 1
    ! h                                  @   
    !          |-----*-----|-----*-----|-----*-----|-----*-----|
    ! z_rho         k-1           k         k+1         k+2 
    ! z_kv                k-1          k          k+1   
    ! Situation 2
    ! h                                        @
    !          |-----*-----|-----*-----|-----*-----|-----*-----|
    ! z_rh0         k-1           k         k+1         k+2 
    ! z_kv                k-1          k          k+1         k+2

    if (hbl .lt. abs(z_rho(1))) then
       k = sum(maxloc(z_q,z_q<= -1.0_idx * hbl))
       If (hbl <= abs(z_rho(k+1))) then
          ! Situation I
          cff_up =  (z_rho(k+1)+hbl) / ((z_q(k+1)-z_q(k)) * (z_rho(k+1)-z_rho(k)));
          cff_dn = (-1.0_idx*hbl-z_rho(k)) / ( (z_q(k)-z_q(k-1)) * (z_rho(k+1)-z_rho(k)))
          KMp = cff_dn * (KM(k)-KM(k-1)) + cff_up * (KM(k+1)-KM(k))
          KTp = cff_dn * (KT(k)-KT(k-1)) + cff_up * (KT(k+1)-KT(k))
       else
          ! Situation II        
          cff_up =  (z_rho(k+2)+hbl) / ((z_q(k+2)-z_q(k+1)) * (z_rho(k+2)-z_rho(k+1)));
          cff_dn = (-1.0_idx*hbl-z_rho(k+1)) / ( (z_q(k+1)-z_q(k)) * (z_rho(k+2)-z_rho(k+1)))
          KMp = cff_dn * (KM(k+1)-KM(k)) + cff_up * (KM(k+2)-KM(k+1))
          KTp = cff_dn * (KT(k+1)-KT(k)) + cff_up * (KT(k+2)-KT(k+1))
       end if
       !cff_up =  (-1.0*hbl - z_q(k))/ (z_q(k+1)-z_q(k))
       !cff_dn =  (z_q(k+1)+hbl) /  (z_q(k+1)-z_q(k))
       !Kmp = (KM(k+1)-KM(k)) / (z_q(k+1)-z_q(k)) ; KTp = (KT(k+1)-KT(k)) / (z_q(k+1)-z_q(k))
       !Kmh = KM(k+1) * cff_up + KM(k) * cff_dn
       !KTh = KT(k+1) * cff_up + KT(k) * cff_dn
       KMh = KM(k)  + KMp*(-1.0*hbl-z_q(k))
       KTh = KT(k)  + KTp*(-1.0*hbl-z_q(k))
       Gm1 = KMh / (hbl * w_m + tiny)
       dGm1ds = min(0.0_idx,-1.0_idx * KMp / (w_m+tiny) - KMh * f1)
       ! Scalar fields 
       !KTh = KT(k) + KTp * (-z_q(k)-hbl);
       Gt1 = KTh / (hbl * w_s + tiny) ;
       dGt1ds =  min(0.0_idx,-1.0_idx * KTp / (w_s+tiny) - KTh * f1)
    else
       write(*,*) "Maximum depth is less than hbl"
    end if
    ! Constants for shape functions as per LMD (17)=================
    a2_m = -2.0_idx + 3.0_idx * Gm1 - dGm1ds;  a3_m = 1.0_idx - 2.0_idx * Gm1 + dGm1ds;
    a2_s = -2.0_idx + 3.0_idx * Gt1 - dGt1ds;  a3_s = 1.0_idx - 2.0_idx * Gt1 + dGt1ds;
 

    L_inv = kappa * Bo_net / (ustar**3+tiny)
    do i = k+1,N
       sl_depth = hbl * epsilon;
       sigma = abs(z_q(i)) / (hbl+tiny);
       ! Velocity scales
       if (L_inv .ge. 0.0_idx) then
          zeta = abs(z_q(i)) * L_inv
       else
          zeta=min(abs(z_q(i)),sl_depth) * L_inv      
       end if
       phim = phi_m(zeta); phis = phi_s(zeta);
       w_m = kappa * ustar / phim;
       w_s = kappa * ustar / phis;
       ! Shape functions LMD (11) 
       Gm = sigma * (1.0+a2_m*sigma+a3_m*sigma**2)
       Gs = sigma * (1.0+a2_s*sigma+a3_s*sigma**2)
       KM(i) = hbl * w_m * Gm;  KT(i) = hbl * w_s * Gs
       KS(i) = hbl * w_s * Gs
    end do
  end subroutine kpp_vdiff
end module kpp_sub

