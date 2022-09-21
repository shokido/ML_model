!********************************************
! Lists of CPP FLAGS                        *
!********************************************
!==============================!
! Atmospheric forcing          !
!==============================!
! Bulk formula
!#undef BULK_FLUXES   ! Use sensible/latent heat flux
#define BULK_FLUXES  ! Use surface air temperature/humidity
#undef BULK_KARA05   ! Use NCEP bulk formula
!#define BULK_KARA05 ! Use Kara et al. (2005) 
!
! Longwave
#undef LONGWAVE_DOWN   ! Use net longwave radiation
!#define LONGWAVE_DOWN ! Use downwelling longwave radiation
!
! Wind speed
!#undef WSPEED_IN       ! Not use wind speed data
#define  WSPEED_IN    ! Use wind speed data
#define WSTRESS_IN
  !
! Diurnal
#undef DIURNAL         ! Not add dirunal cycle
!#define DIURNAL       ! Add dirunal cycle
!
!==============================!
! Mixing scheme                !
!==============================! 			 
#define NNF    ! Use Nakanishi-Niino-Furuichi-Hibiya scheme
!#define KPP   ! Use KPP scheme
!#define KC    ! Use Kantha-Clayson scheme

!==============================!
! Bottom condition
!==============================!
!#undef zero_bottom   
#define zero_bottom   

