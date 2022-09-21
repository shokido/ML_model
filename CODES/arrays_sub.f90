module arrays_sub
  ! Module for 1-D ML arrays
  use param
  implicit none
  !--------------------------------
  ! Input oceanic data
  integer :: ntime_ocn    
  real(idx),allocatable :: lon_grd(:),lat_grd(:),z_rho(:),z_q(:)
  real(idx),allocatable :: cor(:)
  !--------------------------------
  ! 1D oceanic array
  real(idx),allocatable :: temp_next(:),salt_next(:),qq_next(:),l_next(:)
  real(idx),allocatable :: u_next(:),v_next(:)
  real(idx),allocatable :: bvf_1d(:),shear_1d(:)
  real(idx),allocatable :: km_1d(:),kt_1d(:),ks_1d(:),kq_1d(:)
  real(idx),allocatable :: penet(:)
contains

  subroutine deallocate_ocn_1d_arrays(temp,salt,u,v,qq,l)
    implicit none
    real(idx),allocatable :: temp(:),salt(:)
    real(idx),allocatable :: u(:),v(:)
    real(idx),allocatable :: qq(:),l(:)
    deallocate(temp) ; deallocate(salt)
    deallocate(u) ; deallocate(v)
    deallocate(qq) ;  deallocate(l)
  end subroutine deallocate_ocn_1d_arrays
  ! Mixing 1-D array
  subroutine allocate_mix_1d_arrays(nz_q,bvf,shear,km,kt,ks,kq)
    implicit none
    integer,intent(in) :: nz_q ! (nz_q = nz-1)
    real(idx),allocatable :: bvf(:),shear(:),km(:),kt(:),ks(:),kq(:)
    allocate(bvf(nz_q)) ; allocate(shear(nz_q))
    allocate(km(nz_q)) ; allocate(kt(nz_q))
    allocate(ks(nz_q)) ; allocate(kq(nz_q))
    bvf = 0.0_idx ; shear=0.0_idx
    km = 0.0_idx  ; kt = 0.0_idx
    ks = 0.0_idx  ; kq = 0.0_idx
  end subroutine allocate_mix_1d_arrays
  subroutine deallocate_mix_1d_arrays(bvf,shear,km,kt,ks,kq)
    implicit none
    real(idx),allocatable :: bvf(:),shear(:),km(:),kt(:),ks(:),kq(:)
    deallocate(bvf) ; deallocate(shear)
    deallocate(km) ; deallocate(kt)
    deallocate(ks) ; deallocate(kq)
  end subroutine deallocate_mix_1d_arrays
end module arrays_sub
